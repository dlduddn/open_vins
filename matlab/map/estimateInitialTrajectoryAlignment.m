function [x0, info] = estimateInitialTrajectoryAlignment(cfg, mapDB, meas, dSE2, len)
%ESTIMATEINITIALTRAJECTORYALIGNMENT 초기 여러 프레임으로 SE(2) 초기 자세를 추정합니다.
%
% 입력:
%   cfg   - 초기화 및 탐색 설정값입니다.
%   mapDB - 샘플링된 지도 곡선 점입니다. 첫 두 열은 XY 좌표입니다.
%   meas  - 쿼리 궤적의 프레임별 곡선 측정값입니다.
%   dSE2  - 프레임 간 VIO 오도메트리 증가량 [dx; dy; dyaw]입니다.
%   len   - 선택 사항이며, 고려할 프레임 개수입니다.
%
% 출력:
%   x0    - 지도 좌표계에서 추정한 초기 자세 [x; y; yaw]입니다.
%   info  - 시각화, 점수 확인, 디버깅을 위한 진단 정보입니다.
%
% 단일 프레임 곡선 매칭은 직선 도로에서 약합니다. 차선 방향으로의
% 병진 이동은 관측성이 낮기 때문입니다. 이 초기화기는 초기 자세 후보를
% 초기 VIO 오도메트리로 여러 프레임에 전파한 뒤, 사용 가능한 프레임들의
% 곡선-지도 가능도를 누적해서 점수를 계산합니다. 약한 사전값은
% 관측성이 낮은 종방향 성분이 추론 잡음에 과적합되는 것을 막아 줍니다.

    prior = [
        getScalar(cfg, 'initPriorX', 0.0);
        getScalar(cfg, 'initPriorY', 0.0);
        deg2rad(getScalar(cfg, 'initPriorYawDeg', 0.0))
    ];
    x0 = prior;

    % 이후 단계에서 후보 점수를 계산하지 못하더라도 호출자가 사전값을
    % 받을 수 있도록, 유효한 대체 결과로 진단 정보를 초기화합니다.
    info = struct('frameIdx', 1, 'x0', x0, 'prior', prior, ...
        'selectedFrameIdx', [], 'logLikelihood', -inf);

    if nargin < 5 || isempty(len)
        len = numel(meas);
    end

    % 누적 오도메트리 증가량을 첫 프레임 기준 상대 자세로 변환한 뒤,
    % 사용 가능한 베지어 쿼리 곡선이 있는 초기 프레임만 남깁니다.
    relPose = relativeTrajectory(dSE2, len);
    [frameIdx, frames] = selectInitialFrames(cfg, meas, len);
    if isempty(frameIdx)
        warning('estimateInitialTrajectoryAlignment:NoUsableQuery', ...
            'No query Bezier curve survived data association. Returning prior.');
        return;
    end

    % 사전값, 오도메트리, 선택된 쿼리 곡선이 만드는 탐색 영역으로 지도를
    % 잘라냅니다. 이렇게 하면 최근접 곡선 가능도 계산 비용이 줄어듭니다.
    mapCrop = cropMapForTrajectorySearch(cfg, mapDB, frames, frameIdx, relPose, prior);
    if isempty(mapCrop)
        mapCrop = mapDB;
    end
    mapIndex = prepareCurveMapIndex(cfg, mapCrop);

    % 거친 격자 탐색은 사전값을 중심으로 한 전체 불확실성 영역을 탐색합니다.
    searchX = getScalar(cfg, 'initSearchX', 10.0);
    searchY = getScalar(cfg, 'initSearchY', 10.0);
    searchYaw = deg2rad(getScalar(cfg, 'initSearchYawDeg', 30.0));

    coarseX = makeRange(prior(1), searchX, getScalar(cfg, 'initCoarseStepXY', 1.0));
    coarseY = makeRange(prior(2), searchY, getScalar(cfg, 'initCoarseStepXY', 1.0));
    coarseYaw = makeRange(prior(3), searchYaw, deg2rad(getScalar(cfg, 'initCoarseStepYawDeg', 2.0)));
    coarseBest = searchInitialGrid(cfg, mapIndex, frames, frameIdx, relPose, coarseX, coarseY, coarseYaw, prior);

    % 정밀 격자 탐색은 거친 탐색의 최적 자세 주변을 더 작은 간격으로 탐색합니다.
    fineRadiusXY = getScalar(cfg, 'initFineSerachXY', 1.0);
    fineRadiusYaw = deg2rad(getScalar(cfg, 'initFineSearchYawDeg', 2.0));
    fineX = makeRange(coarseBest.x(1), fineRadiusXY, getScalar(cfg, 'initFineStepXY', 0.25));
    fineY = makeRange(coarseBest.x(2), fineRadiusXY, getScalar(cfg, 'initFineStepXY', 0.25));
    fineYaw = makeRange(coarseBest.x(3), fineRadiusYaw, deg2rad(getScalar(cfg, 'initFineStepYawDeg', 0.5)));
    if getLogical(cfg, 'initUseFineSearch', true)
        best = searchInitialGrid(cfg, mapIndex, frames, frameIdx, relPose, fineX, fineY, fineYaw, prior);
    else
        best = coarseBest;
    end

    if ~isfinite(best.logL)
        warning('estimateInitialTrajectoryAlignment:NoFiniteLikelihood', ...
            'All initialization candidates failed map association. Returning prior.');
        x0 = prior;
    else
        x0 = [best.x(1); best.x(2); wrapAngle(best.x(3))];
    end

    % 이후 그리기와 디버깅에 사용할 변환된 쿼리 샘플을 모읍니다.
    sourceXY = collectTrajectorySourceXY(frames, frameIdx, relPose);
    sourceControlXY = collectTrajectoryControlXY(frames, frameIdx, relPose);

    % 초기화기를 다시 실행하지 않고도 사전값, 거친 탐색, 최종 정렬 결과를
    % 시각적으로 비교할 수 있도록 중간 상태를 충분히 저장합니다.
    info.frameIdx = 1;
    info.x0 = x0;
    info.prior = prior;
    info.logLikelihood = best.logL;
    info.coarseX0 = coarseBest.x(:);
    info.coarseLogLikelihood = coarseBest.logL;
    info.selectedFrameIdx = frameIdx;
    info.numSourceCurves = sum(arrayfun(@(f) numel(f.usableIdx), frames));
    info.numSourcePoints = size(sourceXY, 1);
    info.numMapPoints = size(mapCrop, 1);
    info.sourceCurves = [];
    info.sourceXY = sourceXY;
    info.sourceControlXY = sourceControlXY;
    info.mapCrop = mapCrop(:, 1:2);
    info.priorAlignedXY = transformPoints(sourceXY, prior);
    info.coarseAlignedXY = transformPoints(sourceXY, coarseBest.x(:));
    info.alignedXY = transformPoints(sourceXY, x0);

    % 최종 자세에서 선택 프레임들을 다시 평가해 잔차 요약값을 계산합니다.
    evalInfo = evaluateSelectedFrames(cfg, mapIndex, frames, frameIdx, relPose, x0);
    info.likelihoodInfo = evalInfo;
    info.meanResidual = evalInfo.meanResidual;

    fprintf(['Initial trajectory likelihood: frames=%d, idx=[%d..%d], ' ...
             'x0=[%.3f %.3f %.3fdeg], logL=%.2f, curves=%d, samples=%d\n'], ...
        numel(frameIdx), frameIdx(1), frameIdx(end), x0(1), x0(2), ...
        rad2deg(x0(3)), best.logL, info.numSourceCurves, size(sourceXY, 1));
end

function relPose = relativeTrajectory(dSE2, len)
%RELATIVETRAJECTORY 국소 SE(2) 오도메트리를 첫 프레임 기준 자세로 적분합니다.
    relPose = zeros(3, len);
    if isempty(dSE2)
        return;
    end

    for t = 2:len
        % dSE2(:, t-1)는 이전 몸체 좌표계 기준이므로, 누적된 상대 yaw만큼
        % 병진 증가량을 회전해서 첫 프레임 기준 좌표로 누적합니다.
        dx = dSE2(1, t - 1);
        dy = dSE2(2, t - 1);
        dyaw = dSE2(3, t - 1);

        c = cos(relPose(3, t - 1));
        s = sin(relPose(3, t - 1));
        relPose(1, t) = relPose(1, t - 1) + c * dx - s * dy;
        relPose(2, t) = relPose(2, t - 1) + s * dx + c * dy;
        relPose(3, t) = wrapAngle(relPose(3, t - 1) + dyaw);
    end
end

function [selectedIdx, selectedFrames] = selectInitialFrames(cfg, meas, len)
%SELECTINITIALFRAMES 베지어 대응 관계를 통과한 초기 프레임을 선택합니다.
    windowFrames = min([len, numel(meas), round(getScalar(cfg, 'initTrajectoryWindowFrames', 400))]);
    maxFrames = max(1, round(getScalar(cfg, 'initTrajectoryMaxFrames', 12)));

    usableIdx = [];
    frameCache = struct([]);
    for k = 1:windowFrames
        frame = buildQueryBezierFrame(cfg, meas(k));
        if frame.hasUsableCurves
            usableIdx(end + 1) = k; %#ok<AGROW>
            if isempty(frameCache)
                frameCache = frame;
            else
                frameCache(end + 1) = frame; %#ok<AGROW>
            end
        end
    end

    if isempty(usableIdx)
        selectedIdx = [];
        selectedFrames = [];
        return;
    end

    if numel(usableIdx) > maxFrames
        % 초기화 비용을 제한하면서도 궤적 기준 간격을 유지하도록,
        % 사용 가능한 프레임 구간 전체에 고르게 퍼지게 선택합니다.
        pick = unique(round(linspace(1, numel(usableIdx), maxFrames)));
        selectedIdx = usableIdx(pick);
        selectedFrames = frameCache(pick);
    else
        selectedIdx = usableIdx;
        selectedFrames = frameCache;
    end
end

function best = searchInitialGrid(cfg, mapIndex, frames, frameIdx, relPose, xVals, yVals, yawVals, prior)
%SEARCHINITIALGRID 선택된 프레임들에서 모든 초기 자세 후보의 점수를 계산합니다.
    best.x = [0, 0, 0];
    best.logL = -inf;

    [gridX, gridY] = ndgrid(xVals, yVals);
    xyCount = numel(gridX);
    batchSize = max(1, round(getScalar(cfg, 'initBatchSize', 512)));

    for iyaw = 1:numel(yawVals)
        yaw = yawVals(iyaw);
        states0 = [gridX(:).'; gridY(:).'; repmat(yaw, 1, xyCount)];

        for first = 1:batchSize:xyCount
            last = min(first + batchSize - 1, xyCount);
            batchStates0 = states0(:, first:last);

            % 약한 사전값 점수에서 시작한 뒤, 초기 자세를 각 선택 프레임으로
            % 전파해 얻은 지도 가능도를 모두 더합니다.
            score = priorLogScore(cfg, batchStates0, prior);

            for j = 1:numel(frames)
                statesFrame = propagateInitialBatch(batchStates0, relPose(:, frameIdx(j)));
                score = score + curveMapLikelihood(cfg, mapIndex, frames(j), statesFrame);
            end

            [candidateLogL, idx] = max(score);
            if candidateLogL > best.logL
                best.logL = candidateLogL;
                best.x = batchStates0(:, idx).';
            end
        end
    end
end

function score = priorLogScore(cfg, states, prior)
%PRIORLOGSCORE 관측성이 낮은 축을 정규화하는 가우시안 로그 사전값입니다.
    if ~getLogical(cfg, 'initUsePriorScore', true)
        score = zeros(1, size(states, 2));
        return;
    end

    stdX = max(getScalar(cfg, 'initPriorStdX', inf), eps);
    stdY = max(getScalar(cfg, 'initPriorStdY', inf), eps);
    stdYaw = max(deg2rad(getScalar(cfg, 'initPriorStdYawDeg', inf)), eps);

    dx = states(1, :) - prior(1);
    dy = states(2, :) - prior(2);
    dyaw = wrapAngle(states(3, :) - prior(3));
    score = -0.5 * ((dx / stdX) .^ 2 + (dy / stdY) .^ 2 + (dyaw / stdYaw) .^ 2);
end

function statesFrame = propagateInitialBatch(states0, rel)
%PROPAGATEINITIALBATCH 초기 자세 후보를 선택된 프레임의 자세로 전파합니다.
    c = cos(states0(3, :));
    s = sin(states0(3, :));

    statesFrame = states0;
    statesFrame(1, :) = states0(1, :) + c * rel(1) - s * rel(2);
    statesFrame(2, :) = states0(2, :) + s * rel(1) + c * rel(2);
    statesFrame(3, :) = wrapAngle(states0(3, :) + rel(3));
end

function evalInfo = evaluateSelectedFrames(cfg, mapIndex, frames, frameIdx, relPose, x0)
%EVALUATESELECTEDFRAMES 프레임별 매칭 개수와 잔차 요약값을 계산합니다.
    matched = zeros(1, numel(frames));
    residual = inf(1, numel(frames));

    for j = 1:numel(frames)
        stateFrame = propagateInitialBatch(x0(:), relPose(:, frameIdx(j)));
        [~, info] = curveMapLikelihood(cfg, mapIndex, frames(j), stateFrame);
        if isfield(info, 'matchedCurveCount') && ~isempty(info.matchedCurveCount)
            matched(j) = info.matchedCurveCount(1);
        end
        if isfield(info, 'meanResidual') && ~isempty(info.meanResidual)
            residual(j) = info.meanResidual(1);
        end
    end

    evalInfo = struct( ...
        'matchedCurveCount', matched, ...
        'meanResidualByFrame', residual, ...
        'meanResidual', mean(residual(isfinite(residual)), 'omitnan'));
end

function mapCrop = cropMapForTrajectorySearch(cfg, mapDB, frames, frameIdx, relPose, prior)
%CROPMAPFORTRAJECTORYSEARCH 초기화 탐색 영역 안의 지도 점만 남깁니다.
    if isempty(mapDB) || size(mapDB, 2) < 2
        mapCrop = zeros(0, size(mapDB, 2));
        return;
    end

    sourceXY = collectTrajectorySourceXY(frames, frameIdx, relPose);
    if isempty(sourceXY)
        mapCrop = mapDB;
        return;
    end

    searchX = getScalar(cfg, 'initSearchX', 10.0);
    searchY = getScalar(cfg, 'initSearchY', 10.0);
    margin = getScalar(cfg, 'initMapCropMargin', 20.0);

    % 원본 궤적은 아직 첫 프레임 좌표계에 있으므로, 여기서는 사전값의
    % 병진 성분과 설정된 탐색 반경만 사용하면 됩니다.
    xMin = min(sourceXY(:, 1)) + prior(1) - searchX - margin;
    xMax = max(sourceXY(:, 1)) + prior(1) + searchX + margin;
    yMin = min(sourceXY(:, 2)) + prior(2) - searchY - margin;
    yMax = max(sourceXY(:, 2)) + prior(2) + searchY + margin;

    keep = mapDB(:, 1) >= xMin & mapDB(:, 1) <= xMax & ...
           mapDB(:, 2) >= yMin & mapDB(:, 2) <= yMax & ...
           all(isfinite(mapDB(:, 1:2)), 2);
    mapCrop = mapDB(keep, :);
end

function sourceXY = collectTrajectorySourceXY(frames, frameIdx, relPose)
%COLLECTTRAJECTORYSOURCEXY 쿼리 샘플 점을 첫 프레임 좌표계로 모읍니다.
    sourceXY = zeros(0, 2);
    for j = 1:numel(frames)
        rel = relPose(:, frameIdx(j));
        for i = frames(j).usableIdx(:).'
            sourceXY = [sourceXY; transformBodyToInitial(frames(j).curves(i).sampleXY, rel)]; %#ok<AGROW>
        end
    end
end

function sourceXY = collectTrajectoryControlXY(frames, frameIdx, relPose)
%COLLECTTRAJECTORYCONTROLXY 베지어 제어점을 NaN 구분자와 함께 모읍니다.
    sourceXY = zeros(0, 2);
    for j = 1:numel(frames)
        rel = relPose(:, frameIdx(j));
        for i = frames(j).usableIdx(:).'
            sourceXY = [sourceXY; transformBodyToInitial(frames(j).curves(i).controlXY, rel); NaN, NaN]; %#ok<AGROW>
        end
    end
end

function outXY = transformBodyToInitial(bodyXY, rel)
%TRANSFORMBODYTOINITIAL 몸체 좌표계의 곡선 점을 첫 프레임 좌표계로 옮깁니다.
    c = cos(rel(3));
    s = sin(rel(3));
    outXY = bodyXY * [c, s; -s, c] + rel(1:2).';
end

function outXY = transformPoints(sourceXY, x)
%TRANSFORMPOINTS 첫 프레임 기준 궤적 점에 SE(2) 자세를 적용합니다.
    c = cos(x(3));
    s = sin(x(3));
    outXY = sourceXY * [c, s; -s, c] + x(1:2).';
end

function values = makeRange(center, radius, step)
%MAKERANGE 중심값 주변에 대칭 격자를 만듭니다.
    values = center + (-radius:step:radius);
    if isempty(values)
        values = center;
    end
end
