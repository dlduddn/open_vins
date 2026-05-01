function [logL, info] = curveMapLikelihood(cfg, mapDB, queryInput, xStates)
%CURVEMAPLIKELIHOOD Pose likelihood from Query Bezier curves to map centerlines.
%
% xStates is [x; y; yaw] in map-local coordinates. Each query curve is first
% represented as a cubic Bezier centerline and sampled along the curve; raw
% control points are not used as matching points.
%
% 이 함수는 각 후보 pose가 현재 query Bezier curve 관측을 지도 중심선에
% 얼마나 잘 정렬시키는지 log-likelihood로 평가한다. 값이 0에 가까울수록
% 좋은 정렬이고, 큰 음수일수록 관측과 지도가 잘 맞지 않는 pose다.

    % 입력 pose를 내부 표준 형태인 3xN 행렬로 맞춘다. 단일 pose vector는
    % column vector로 만들고, Nx3 입력은 3xN으로 transpose한다.
    if isvector(xStates)
        xStates = xStates(:);
    end
    if size(xStates, 1) ~= 3 && size(xStates, 2) == 3
        xStates = xStates.';
    end

    nState = size(xStates, 2);

    % 기본 log-likelihood는 0이다. 관측이 없거나 usable curve가 없어서
    % 바로 return되면 모든 pose가 같은 점수를 받으므로 PF weight는 변하지 않는다.
    logL = zeros(1, nState);

    % queryInput은 이미 buildQueryBezierFrame()을 거친 frame일 수도 있고,
    % raw measurement일 수도 있다. raw이면 여기서 cubic Bezier frame으로 변환한다.
    frame = queryInput;
    if ~isfield(frame, 'curves')
        frame = buildQueryBezierFrame(cfg, queryInput);
    end

    % info는 likelihood 계산 결과를 해석하기 위한 진단용 구조체다.
    % 각 pose가 gate를 통과했는지, 몇 개 curve가 matched 되었는지 등을 담는다.
    info = emptyInfo(nState);
    info.hasMeasurement = frame.hasUsableCurves;
    info.usableCurveCount = numel(frame.usableIdx);
    info.frame = frame;

    % usable curve가 없으면 이 frame은 pose를 구분하는 정보를 주지 않는다고 보고
    % logL=0을 유지한 채 반환한다.
    if ~frame.hasUsableCurves
        return;
    end

    % mapDB를 nearest-neighbor 검색에 적합한 index 구조로 맞춘다.
    % 이미 index가 들어온 경우에는 재생성하지 않고 그대로 사용한다.
    mapIndex = asMapIndex(cfg, mapDB);
    if isempty(mapIndex.xy)
        % 지도 중심선 점이 없으면 관측과 비교할 대상이 없으므로 모든 pose에
        % miss penalty를 부여한다.
        logL(:) = getScalar(cfg, 'likelihoodMissLogPenalty', -60.0);
        info.hasMeasurement = false;
        return;
    end

    % association과 likelihood를 제어하는 설정값들이다.
    % sigma는 거리 residual의 scale, maxDist는 matched sample 인정 거리,
    % minMatchFraction은 한 curve가 matched로 인정되기 위한 sample 비율이다.
    sigma = max(getScalar(cfg, 'likelihoodSigma', 0.75), eps);
    robustScale = max(getScalar(cfg, 'likelihoodRobustScale', 1.0), eps);
    maxDist = getScalar(cfg, 'assocMaxCurveMapDist', 2.5);
    maxD2 = maxDist^2;
    minMatchFraction = getScalar(cfg, 'assocMinCurveMatchFraction', 0.45);
    minMatchedCurves = getScalar(cfg, 'assocMinMatchedCurves', 1);
    maxEffSamples = getScalar(cfg, 'likelihoodMaxEffectiveSamples', 25);
    curveMissPenalty = -abs(getScalar(cfg, 'likelihoodCurveMissLogPenalty', 30.0));
    maxCurvePenalty = abs(getScalar(cfg, 'likelihoodMaxCurvePenalty', abs(curveMissPenalty)));
    inlierReward = getScalar(cfg, 'likelihoodInlierReward', 0.0);
    missPenalty = -abs(getScalar(cfg, 'likelihoodMissLogPenalty', 60.0));
    useRobustLikelihood = getLogical(cfg, 'likelihoodUseRobust', true);
    useBodyGate = getLogical(cfg, 'assocUseBodyRoadGate', true);
    useCurveCorrespondence = getLogical(cfg, 'assocUseCurveMapCorrespondenceRule', true);
    fovGate = computeCameraFovGate(cfg);

    % frame 내 모든 curve가 아니라 전처리 규칙을 통과한 usable curve만 사용한다.
    curves = frame.curves(frame.usableIdx);

    % NaN/Inf가 포함된 pose는 geometry 계산이 불가능하므로 miss 처리한다.
    finiteState = all(isfinite(xStates), 1);
    logL(~finiteState) = missPenalty;

    if useBodyGate
        % body-road gate는 pose 자체가 지도 도로 주변에 있는지 먼저 확인한다.
        % 명백히 도로 밖에 있는 pose는 비싼 curve matching 전에 제거한다.
        [roadOk, lateral, lateralLimit] = bodyRoadGateBatch(cfg, mapIndex, xStates);
        info.roadGateOk = roadOk;
        info.bodyLateralDistance = lateral;
        info.bodyLateralLimit = lateralLimit;
    else
        % gate를 끄면 finite pose는 모두 curve matching 평가 대상으로 둔다.
        roadOk = true(1, nState);
        info.roadGateOk = roadOk;
    end

    % 실제 curve-map likelihood를 계산할 pose만 선택한다.
    evalMask = finiteState & roadOk;
    logL(~evalMask) = missPenalty;
    if ~any(evalMask)
        return;
    end

    evalIdx = find(evalMask);
    if fovGate.enabled && strcmp(fovGate.mode, "candidate")
        % candidate FOV mode:
        % 각 pose 주변에서 카메라 FOV 안에 들어올 수 있는 map point 후보를 먼저
        % 줄인 뒤, 그 후보 집합 안에서 거리 기반 likelihood를 계산한다.
        [stateLogL, matchedCount, meanResidual] = evaluateStatesInFov( ...
            curves, mapIndex, xStates(:, evalIdx), sigma, maxD2, maxDist, ...
            minMatchFraction, maxEffSamples, fovGate, robustScale, ...
            curveMissPenalty, maxCurvePenalty, inlierReward, useRobustLikelihood, ...
            useCurveCorrespondence);
    else
        % 기본 mode:
        % query sample을 pose별 map 좌표로 변환한 뒤 전체 map index에서 최근접점을
        % 찾고, 필요하면 nearest point가 FOV 안에 있는지 추가로 검사한다.
        [stateLogL, matchedCount, meanResidual] = evaluateStates( ...
            curves, mapIndex, xStates(:, evalIdx), sigma, maxD2, maxDist, ...
            minMatchFraction, maxEffSamples, fovGate, robustScale, ...
            curveMissPenalty, maxCurvePenalty, inlierReward, useRobustLikelihood, ...
            useCurveCorrespondence);
    end

    % 디버깅/로그용 결과를 원래 pose index 위치에 되돌려 기록한다.
    info.matchedCurveCount(evalIdx) = matchedCount;
    info.meanResidual(evalIdx) = meanResidual;

    % 최소 matched curve 수를 만족한 pose만 association 성공으로 인정한다.
    % 실패한 pose는 missPenalty, 성공한 pose는 거리 기반 stateLogL을 사용한다.
    if useCurveCorrespondence
        associated = matchedCount >= minMatchedCurves;
    else
        associated = true(size(matchedCount));
    end
    if useRobustLikelihood
        frameMissPenalty = min(missPenalty, curveMissPenalty * numel(curves));
    else
        frameMissPenalty = missPenalty;
    end
    logL(evalIdx(~associated)) = frameMissPenalty;
    logL(evalIdx(associated)) = stateLogL(associated);
end

function mapIndex = asMapIndex(cfg, mapInput)
    % mapInput이 이미 prepareCurveMapIndex() 결과라면 재사용하고, raw mapDB이면
    % 여기서 index를 만든다.
    if isstruct(mapInput) && isfield(mapInput, 'xy') && isfield(mapInput, 'width')
        mapIndex = mapInput;
    else
        mapIndex = prepareCurveMapIndex(cfg, mapInput);
    end
end

function [stateLogL, matchedCount, meanResidual] = evaluateStatesInFov( ...
        curves, mapIndex, xStates, sigma, maxD2, maxDist, minMatchFraction, ...
        maxEffSamples, fovGate, robustScale, curveMissPenalty, maxCurvePenalty, ...
        inlierReward, useRobustLikelihood, useCurveCorrespondence)

    nState = size(xStates, 2);
    nCurve = numel(curves);

    % pose별 누적 likelihood와 진단값이다. 이 함수는 각 pose를 독립적으로
    % 순회하면서 FOV 안의 map 후보만 대상으로 association을 수행한다.
    if useRobustLikelihood
        stateLogL = curveMissPenalty * nCurve * ones(1, nState);
    else
        stateLogL = zeros(1, nState);
    end
    matchedCount = zeros(1, nState);
    residualSum = zeros(1, nState);
    residualCount = zeros(1, nState);

    % FOV 최대 range에 association margin(maxDist)을 더해 검색 반경을 잡는다.
    % 경계 근처의 map point가 미리 잘려나가지 않게 하기 위한 여유다.
    searchRadius = fovGate.maxRange + maxDist;
    if mapIndex.hasKdTree && exist('rangesearch', 'file') == 2
        candidateCells = rangesearch(mapIndex.searcher, xStates(1:2, :).', searchRadius);
    else
        % rangesearch를 사용할 수 없으면 모든 map point를 후보로 사용한다.
        % 정확도는 유지되지만 계산량은 커진다.
        allIdx = 1:size(mapIndex.xy, 1);
        candidateCells = repmat({allIdx}, 1, nState);
    end

    for j = 1:nState
        idx = candidateCells{j};
        if isempty(idx)
            continue;
        end

        % 현재 pose 주변 map 후보를 body 좌표계로 변환해서 카메라 FOV 안에
        % 들어오는 후보만 남긴다.
        c = cos(xStates(3, j));
        s = sin(xStates(3, j));
        rel = mapIndex.xy(idx, :) - xStates(1:2, j).';
        mapBody = [rel(:, 1) * c + rel(:, 2) * s, ...
                  -rel(:, 1) * s + rel(:, 2) * c];
        fovKeep = pointsInsideCameraFov(mapBody, fovGate, maxDist);
        if ~any(fovKeep)
            continue;
        end

        candidateMapXY = mapIndex.xy(idx(fovKeep), :);
        R = [c, -s; s, c];
        t = xStates(1:2, j).';

        for i = 1:numel(curves)
            sampleXY = curves(i).sampleXY;

            % query curve sample 중 카메라 FOV 안에 실제로 보이는 sample만 평가한다.
            sampleKeep = pointsInsideCameraFov(sampleXY, fovGate, 0.0);
            if ~any(sampleKeep)
                continue;
            end

            % 보이는 query sample을 현재 pose로 map 좌표계에 투영한 뒤, 미리 줄여둔
            % FOV map 후보 집합 안에서 최근접 거리를 구한다.
            visibleSampleXY = sampleXY(sampleKeep, :);
            queryXY = visibleSampleXY * R' + t;
            d2 = nearestSubsetSquaredDistances(queryXY, candidateMapXY);
            dist = sqrt(d2);

            % 한 curve의 sample 중 충분한 비율이 maxDist 안에 들어와야 matched로 본다.
            matchFraction = mean(dist <= sqrt(maxD2));
            if useCurveCorrespondence && matchFraction < minMatchFraction
                continue;
            end

            % 거리 residual을 maxDist^2로 clipping해 outlier 하나가 likelihood를
            % 무한히 망가뜨리지 않게 한다.
            d2Clipped = min(d2, maxD2);
            nEff = min(numel(d2Clipped), maxEffSamples);
            curveLogL = computeCurveLogLikelihood(d2Clipped, matchFraction, ...
                sigma, robustScale, nEff, maxCurvePenalty, inlierReward, useRobustLikelihood);
            if useRobustLikelihood
                stateLogL(j) = stateLogL(j) - curveMissPenalty + curveLogL;
            else
                stateLogL(j) = stateLogL(j) + curveLogL;
            end
            if matchFraction >= minMatchFraction || ~useCurveCorrespondence
                matchedCount(j) = matchedCount(j) + 1;
            end

            % meanResidual은 score 계산용이 아니라 로그/디버깅용 진단값이다.
            residualSum(j) = residualSum(j) + sum(dist);
            residualCount(j) = residualCount(j) + numel(dist);
        end
    end

    meanResidual = inf(1, nState);
    hasResidual = residualCount > 0;
    meanResidual(hasResidual) = residualSum(hasResidual) ./ residualCount(hasResidual);
end

function [stateLogL, matchedCount, meanResidual] = evaluateStates( ...
        curves, mapIndex, xStates, sigma, maxD2, maxDist, minMatchFraction, ...
        maxEffSamples, fovGate, robustScale, curveMissPenalty, maxCurvePenalty, ...
        inlierReward, useRobustLikelihood, useCurveCorrespondence)

    nState = size(xStates, 2);
    nCurve = numel(curves);

    % pose별 누적 likelihood와 진단값. 여러 curve의 score가 stateLogL에 더해진다.
    if useRobustLikelihood
        stateLogL = curveMissPenalty * nCurve * ones(1, nState);
    else
        stateLogL = zeros(1, nState);
    end
    matchedCount = zeros(1, nState);
    residualSum = zeros(1, nState);
    residualCount = zeros(1, nState);

    c = cos(xStates(3, :));
    s = sin(xStates(3, :));
    for i = 1:numel(curves)
        sampleXY = curves(i).sampleXY;

        % body 좌표계의 Bezier sample point를 각 pose 후보의 map 좌표로 변환한다.
        % sampleXY는 nSamplex2, c/s/xStates는 1xN이므로 implicit expansion으로
        % nSamplexN 크기의 qx/qy가 만들어진다.
        qx = sampleXY(:, 1) .* c - sampleXY(:, 2) .* s + xStates(1, :);
        qy = sampleXY(:, 1) .* s + sampleXY(:, 2) .* c + xStates(2, :);
        nSample = size(sampleXY, 1);

        % 모든 pose의 모든 transformed sample을 한 번에 펼쳐 최근접 map point를 찾는다.
        % 반환된 벡터는 아래에서 다시 nSample x nState로 복원한다.
        [d2, nnIdx] = nearestMapSquaredDistances(mapIndex, [qx(:), qy(:)]);
        if fovGate.enabled
            % FOV gate가 켜져 있으면 nearest map point가 해당 pose의 카메라 시야
            % 안에 있을 때만 association을 허용한다.
            fovKeep = nearestPointsInsideFov(mapIndex, xStates, nnIdx, nSample, fovGate, maxDist);
            d2(~fovKeep) = inf;
        end
        d2 = reshape(d2, size(sampleXY, 1), nState);

        dist = sqrt(d2);

        % 각 pose에서 이 curve의 sample 중 maxDist 이내에 들어온 비율을 계산한다.
        matchFraction = mean(dist <= maxDist, 1);
        matched = matchFraction >= minMatchFraction;
        scored = matched | ~useCurveCorrespondence;

        if ~any(scored)
            continue;
        end

        % scored pose에 대해서만 likelihood를 누적한다. residual은 maxD2로 clipping하고,
        % sample 수 효과는 maxEffSamples로 제한한다.
        d2Clipped = min(d2, maxD2);
        nEff = min(size(d2Clipped, 1), maxEffSamples);
        curveLogL = computeCurveLogLikelihood(d2Clipped, matchFraction, ...
            sigma, robustScale, nEff, maxCurvePenalty, inlierReward, useRobustLikelihood);

        if useRobustLikelihood
            stateLogL(scored) = stateLogL(scored) - curveMissPenalty + curveLogL(scored);
        else
            stateLogL(scored) = stateLogL(scored) + curveLogL(scored);
        end
        matchedCount(matched) = matchedCount(matched) + 1;

        % 평균 residual은 scored pose에 대해서만 누적한다. FOV 밖 point 때문에
        % dist에 inf가 있으면 meanResidual도 inf가 될 수 있다.
        residualSum(scored) = residualSum(scored) + sum(dist(:, scored), 1);
        residualCount(scored) = residualCount(scored) + size(dist, 1);
    end

    meanResidual = inf(1, nState);
    hasResidual = residualCount > 0;
    meanResidual(hasResidual) = residualSum(hasResidual) ./ residualCount(hasResidual);
end

function curveLogL = computeCurveLogLikelihood(d2Clipped, matchFraction, ...
        sigma, robustScale, nEff, maxCurvePenalty, inlierReward, useRobustLikelihood)
    if ~useRobustLikelihood
        normalizedCost = mean(d2Clipped, 1) / (sigma^2);
        curveLogL = -0.5 * nEff * normalizedCost;
        return;
    end

    % Robust Cauchy-style cost. The cap keeps an accepted geometric match
    % better than a complete miss, while still ranking close matches highest.
    scale2 = max((sigma * robustScale)^2, eps);
    rho = log1p(d2Clipped ./ scale2);
    curveLogL = -0.5 * nEff * mean(rho, 1) + inlierReward .* matchFraction;
    curveLogL(~isfinite(curveLogL)) = -abs(maxCurvePenalty);
    curveLogL = max(curveLogL, -abs(maxCurvePenalty));
end

function keep = pointsInsideCameraFov(bodyXY, gate, rangeMargin)
    % body 좌표계 점들이 카메라 FOV 안에 있는지 검사한다.
    % rangeMargin은 association 경계 근처 후보를 조금 더 허용하기 위한 여유 거리다.
    rangeMargin = max(rangeMargin, 0.0);
    relAngle = wrapAngle(atan2(bodyXY(:, 2), bodyXY(:, 1)) - gate.yawOffset);
    range = sqrt(sum(bodyXY .^ 2, 2));

    % 카메라 yaw offset 방향으로의 전방 성분이다. 단순히 body x>0을 보는 대신
    % 카메라가 틀어진 방향을 기준으로 앞쪽 여부를 판단한다.
    forward = bodyXY(:, 1) * cos(gate.yawOffset) + bodyXY(:, 2) * sin(gate.yawOffset);

    keep = forward > 0 & ...
        range >= max(0.0, gate.minRange - rangeMargin) & ...
        range <= gate.maxRange + rangeMargin & ...
        relAngle >= -gate.rightAngle - gate.margin & ...
        relAngle <= gate.leftAngle + gate.margin;
end

function d2 = nearestSubsetSquaredDistances(queryXY, mapXY)
    % 주어진 mapXY 부분집합 안에서 각 query point의 최근접 제곱거리를 계산한다.
    % candidate FOV mode에서 pose별 map 후보가 이미 줄어든 경우 사용한다.
    if isempty(mapXY)
        d2 = inf(size(queryXY, 1), 1);
        return;
    end

    d2 = inf(size(queryXY, 1), 1);
    mapX = mapXY(:, 1).';
    mapY = mapXY(:, 2).';
    chunkSize = 256;
    for first = 1:chunkSize:size(queryXY, 1)
        % 큰 거리 행렬을 한 번에 만들지 않도록 query point를 chunk로 나눠 처리한다.
        last = min(first + chunkSize - 1, size(queryXY, 1));
        dx = queryXY(first:last, 1) - mapX;
        dy = queryXY(first:last, 2) - mapY;
        d2(first:last) = min(dx .* dx + dy .* dy, [], 2);
    end
end

function [ok, lateral, lateralLimit] = bodyRoadGateBatch(cfg, mapIndex, xStates)
    % pose 주변의 도로 중심선을 찾아 차량 body가 도로 폭 안에 있는지 검사한다.
    % 통과하지 못한 pose는 curve association 전에 missPenalty 처리된다.
    nState = size(xStates, 2);
    ok = false(1, nState);
    lateral = inf(1, nState);
    lateralLimit = NaN(1, nState);

    % body 기준 검사 영역은 뒤쪽 back부터 앞쪽 lookahead까지다.
    back = getScalar(cfg, 'assocBodyGateBack', 2.0);
    lookahead = getScalar(cfg, 'assocBodyGateLookahead', 30.0);
    defaultWidth = getScalar(cfg, 'assocDefaultRoadWidth', 6.0);
    widthScale = getScalar(cfg, 'assocRoadWidthScale', 1.0);

    % 검색 반경은 전방 길이와 최대 도로 폭을 함께 고려해 잡는다.
    maxWidth = maxFinite([mapIndex.maxRoadWidth; defaultWidth], defaultWidth);
    searchRadius = hypot(max(back, lookahead), widthScale * maxWidth);

    finiteState = all(isfinite(xStates), 1);
    evalIdx = find(finiteState);
    candidateCells = cell(1, nState);
    if isempty(evalIdx)
        return;
    end
    if mapIndex.hasKdTree && exist('rangesearch', 'file') == 2
        % finite pose에 대해서만 주변 map point 후보를 검색한다.
        evalCandidates = rangesearch(mapIndex.searcher, xStates(1:2, evalIdx).', searchRadius);
        for k = 1:numel(evalIdx)
            candidateCells{evalIdx(k)} = evalCandidates{k};
        end
    else
        % KD-tree를 사용할 수 없으면 모든 map point를 후보로 둔다.
        allIdx = 1:size(mapIndex.xy, 1);
        candidateCells(evalIdx) = {allIdx};
    end

    for j = 1:nState
        idx = candidateCells{j};
        if isempty(idx) || any(~isfinite(xStates(:, j)))
            continue;
        end

        c = cos(xStates(3, j));
        s = sin(xStates(3, j));

        % 후보 map point를 pose 기준 body 좌표계로 변환한다.
        rel = mapIndex.xy(idx, :) - xStates(1:2, j).';
        forward = rel(:, 1) * c + rel(:, 2) * s;
        lateralSigned = -rel(:, 1) * s + rel(:, 2) * c;

        % 차량 근처 앞/뒤 구간 안의 map point만 도로 후보로 본다.
        keep = forward >= -back & forward <= lookahead;
        if ~any(keep)
            continue;
        end

        % 가장 가까운 lateral centerline distance를 도로 중심과의 거리로 사용한다.
        keptIdx = idx(keep);
        lateralAbs = abs(lateralSigned(keep));
        [lateral(j), localIdx] = min(lateralAbs);
        mapIdx = keptIdx(localIdx);

        % 지도 폭 정보가 없거나 유효하지 않으면 기본 도로 폭을 사용한다.
        roadWidth = mapIndex.width(mapIdx);
        if ~isfinite(roadWidth) || roadWidth <= 0
            roadWidth = defaultWidth;
        end

        % lateral distance가 허용 폭 안이면 road gate 통과다.
        lateralLimit(j) = widthScale * roadWidth;
        ok(j) = lateral(j) <= lateralLimit(j);
    end
end

function keep = nearestPointsInsideFov(mapIndex, xStates, nnIdx, nSample, fovGate, rangeMargin)
    % evaluateStates()에서 펼쳐진 query sample들의 nearest map point가 각 pose의
    % 카메라 FOV 안에 있는지 검사한다.
    nState = size(xStates, 2);

    % qx(:), qy(:)가 sample-major 순서로 펼쳐졌으므로 각 row가 어떤 state에
    % 속하는지 같은 순서로 복원한다.
    stateIdx = repelem(1:nState, nSample).';
    nearestXY = mapIndex.xy(nnIdx, :);
    stateXY = xStates(1:2, stateIdx).';
    yaw = xStates(3, stateIdx).';

    % nearest map point를 해당 state의 body 좌표계로 변환한 뒤 FOV 검사를 수행한다.
    rel = nearestXY - stateXY;
    c = cos(yaw);
    s = sin(yaw);
    bodyXY = [rel(:, 1) .* c + rel(:, 2) .* s, ...
             -rel(:, 1) .* s + rel(:, 2) .* c];
    keep = pointsInsideCameraFov(bodyXY, fovGate, rangeMargin);
end

function [d2, idx] = nearestMapSquaredDistances(mapIndex, queryXY)
    % 각 query point에 대해 지도 중심선에서 가장 가까운 point와 제곱거리를 찾는다.
    if mapIndex.hasKdTree
        % KD-tree가 준비되어 있으면 빠른 최근접 검색을 사용한다.
        [idx, dist] = knnsearch(mapIndex.searcher, queryXY, 'K', 1);
        d2 = dist .^ 2;
        return;
    end

    % KD-tree가 없을 때의 fallback. 정확한 brute-force 검색이지만 chunk로 나눠
    % 메모리 사용량을 제한한다.
    d2 = inf(size(queryXY, 1), 1);
    idx = zeros(size(queryXY, 1), 1);
    mapX = mapIndex.xy(:, 1).';
    mapY = mapIndex.xy(:, 2).';
    chunkSize = 512;
    for first = 1:chunkSize:size(queryXY, 1)
        last = min(first + chunkSize - 1, size(queryXY, 1));
        dx = queryXY(first:last, 1) - mapX;
        dy = queryXY(first:last, 2) - mapY;
        [d2(first:last), localIdx] = min(dx .* dx + dy .* dy, [], 2);
        idx(first:last) = localIdx;
    end
end

function info = emptyInfo(nState)
    % likelihood 계산 과정에서 얻은 진단 정보를 담는 구조체의 기본값.
    info = struct( ...
        'hasMeasurement', false, ...
        'usableCurveCount', 0, ...
        'matchedCurveCount', zeros(1, nState), ...
        'meanResidual', inf(1, nState), ...
        'roadGateOk', false(1, nState), ...
        'bodyLateralDistance', inf(1, nState), ...
        'bodyLateralLimit', NaN(1, nState), ...
        'frame', []);
end
