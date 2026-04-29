function [x0, info] = estimateInitialAlignment(cfg, mapDB, meas)
%ESTIMATEINITIALALIGNMENT 간단한 SE(2) 곡선-지도 초기 정렬.
%
% Bezier 측정값을 네 개의 독립적인 점이 아니라 3차 중심선 곡선으로
% 처리한다. 각 곡선을 촘촘히 샘플링한 뒤 지도 중심선 점군과 곡선
% 단위로 비교해 점수를 계산한다:
%
%   p_map = R(x0(3)) * p_body + x0(1:2)

    % 입력 검증 단계에서 실패하면 기본 자세를 반환한다.
    info = struct();
    x0 = [0; 0; 0];

    % 유효한 지도 XY 점만 남긴다. 정렬 비용은 이 중심선 점군을 기준으로
    % 계산된다.
    mapXY = mapDB(:, 1:2);
    mapXY = mapXY(all(isfinite(mapXY), 2), :);
    if isempty(mapXY)
        warning('estimateInitialAlignment:EmptyMap', ...
            'Map DB has no valid XY points. Returning x0=[0;0;0].');
        return;
    end

    % 하나 이상의 사용 가능한 중심선 곡선으로 변환할 수 있는 첫 측정
    % 프레임을 사용한다.
    [frameIdx, sourceCurves] = firstValidCurveMeasurement(cfg, meas);
    if isempty(sourceCurves)
        warning('estimateInitialAlignment:EmptyMeasurement', ...
            'No valid Bezier centerline curve found. Returning x0=[0;0;0].');
        return;
    end

    % 샘플링된 점은 최근접 이웃 정렬 비용 계산에 사용된다. 원래 제어점은
    % 진단과 시각화를 위해 별도로 보관한다.
    sourceXY = collectCurveSamples(sourceCurves);
    sourceControlXY = collectCurveControls(sourceCurves);
    if size(sourceXY, 1) < 2
        warning('estimateInitialAlignment:TooFewPoints', ...
            'Too few sampled Bezier curve points for alignment. Returning x0=[0;0;0].');
        return;
    end

    % 사전 자세(prior)는 격자 탐색의 중심이다. 병진 이동은 미터 단위이고
    % 요각(yaw)은 내부적으로 라디안 단위로 저장한다.
    prior = [
        getScalar(cfg, 'alignPriorX', 0);
        getScalar(cfg, 'alignPriorY', 0);
        deg2rad(getScalar(cfg, 'alignPriorYawDeg', 0))
    ];

    % 사전 자세를 중심으로 하는 탐색 반경이다.
    searchX = getScalar(cfg, 'alignSearchX', 10.0);
    searchY = getScalar(cfg, 'alignSearchY', 10.0);
    searchYaw = deg2rad(getScalar(cfg, 'alignSearchYawDeg', 15.0));

    % 완전 탐색 중 최근접 이웃 계산량을 줄이기 위해, 가능한 변환 범위를
    % 포함하도록 지도를 보수적으로 잘라낸다.
    cropMargin = getScalar(cfg, 'alignMapCropMargin', 20.0);
    mapCrop = cropMap(mapXY, sourceXY, prior, searchX, searchY, cropMargin);
    if isempty(mapCrop)
        mapCrop = mapXY;
    end

    % 큰 최근접 이웃 오차를 잘라내어 일부 나쁜 구간이 점수를 지배하지
    % 않게 한다. 곡선 단위 절사로 가장 나쁜 곡선들도 완화한다.
    maxMatchDist = getScalar(cfg, 'alignMaxMatchDist', 3.0);
    d2Clip = maxMatchDist^2;
    trimFraction = getScalar(cfg, 'alignCurveTrimFraction', 0.8);

    % 거친 탐색은 설정된 전체 탐색 영역을 훑는다.
    coarseStepXY = getScalar(cfg, 'alignCoarseStepXY', 1.0);
    coarseStepYaw = deg2rad(getScalar(cfg, 'alignCoarseStepYawDeg', 2.0));

    coarseX = makeRange(prior(1), searchX, coarseStepXY);
    coarseY = makeRange(prior(2), searchY, coarseStepXY);
    coarseYaw = makeRange(prior(3), searchYaw, coarseStepYaw);
    best = gridSearch(sourceCurves, mapCrop, coarseX, coarseY, coarseYaw, d2Clip, trimFraction);
    coarseBest = best;

    % 정밀 탐색은 거친 탐색에서 찾은 최적값 주변만 더 촘촘히 탐색한다.
    fineRadiusXY = getScalar(cfg, 'alignFineRadiusXY', 1.0);
    fineRadiusYaw = deg2rad(getScalar(cfg, 'alignFineRadiusYawDeg', 2.0));
    fineStepXY = getScalar(cfg, 'alignFineStepXY', 0.25);
    fineStepYaw = deg2rad(getScalar(cfg, 'alignFineStepYawDeg', 0.5));

    fineX = makeRange(best.x(1), fineRadiusXY, fineStepXY);
    fineY = makeRange(best.x(2), fineRadiusXY, fineStepXY);
    fineYaw = makeRange(best.x(3), fineRadiusYaw, fineStepYaw);
    best = gridSearch(sourceCurves, mapCrop, fineX, fineY, fineYaw, d2Clip, trimFraction);

    x0 = [best.x(1); best.x(2); wrapAngle(best.x(3))];
    alignedXY = transformPoints(sourceXY, x0);
    d2 = nearestSquaredDistances(alignedXY, mapCrop);

    % 호출자가 정렬 품질을 시각화하거나 디버깅할 수 있도록 최종 결과와
    % 중간 산출물을 함께 저장한다.
    info.frameIdx = frameIdx;
    info.x0 = x0;
    info.prior = prior;
    info.cost = best.cost;
    info.curveCosts = best.curveCosts;
    info.coarseX0 = coarseBest.x(:);
    info.coarseCost = coarseBest.cost;
    info.coarseCurveCosts = coarseBest.curveCosts;
    info.medianError = sqrt(median(d2));
    info.meanClippedError = sqrt(mean(min(d2, d2Clip)));
    info.numSourceCurves = numel(sourceCurves);
    info.numSourcePoints = size(sourceXY, 1);
    info.numMapPoints = size(mapCrop, 1);
    info.maxMatchDist = maxMatchDist;
    info.sourceCurves = sourceCurves;
    info.sourceXY = sourceXY;
    info.sourceControlXY = sourceControlXY;
    info.mapCrop = mapCrop;
    info.priorAlignedXY = transformPoints(sourceXY, prior);
    info.coarseAlignedXY = transformPoints(sourceXY, coarseBest.x(:));
    info.alignedXY = alignedXY;

    fprintf(['Initial curve alignment: frame=%d, x0=[%.3f %.3f %.3fdeg], ' ...
             'medianNN=%.3fm, curves=%d, samples=%d, map=%d\n'], ...
        info.frameIdx, x0(1), x0(2), rad2deg(x0(3)), ...
        info.medianError, info.numSourceCurves, info.numSourcePoints, info.numMapPoints);
end

function [frameIdx, curves] = firstValidCurveMeasurement(cfg, meas)
    % 측정값을 시간 순서로 훑고, 최소 하나의 곡선을 만들 수 있을 만큼
    % 차체 좌표계 점이 충분한 첫 프레임에서 멈춘다.
    frameIdx = [];
    curves = [];

    for k = 1:numel(meas)
        if isfield(meas(k), 'pointsBody') && size(meas(k).pointsBody, 1) >= 4
            curves = measurementToCurves(cfg, meas(k).pointsBody);
            if ~isempty(curves)
                frameIdx = k;
                return;
            end
        end
    end
end

function curves = measurementToCurves(cfg, pointsBody)
    % 하나의 측정 프레임을 곡선 구조체로 변환한다. 4열이 있으면 그 값을
    % 곡선 ID로 사용해 제어점을 그룹화한다.
    curves = [];
    pointsBody = pointsBody(all(isfinite(pointsBody(:, 1:2)), 2), :);
    if size(pointsBody, 1) < 4
        return;
    end

    if size(pointsBody, 2) >= 4
        curveIds = unique(pointsBody(:, 4), 'stable');
    else
        curveIds = 1;
        pointsBody(:, 4) = 1;
    end

    for i = 1:numel(curveIds)
        id = curveIds(i);
        controlXY = pointsBody(pointsBody(:, 4) == id, 1:2);
        if size(controlXY, 1) < 4
            continue;
        end

        newCurves = controlPointsToCurves(cfg, id, controlXY);
        curves = appendCurves(curves, newCurves);
    end
end

function curves = controlPointsToCurves(cfg, id, controlXY)
    % 각 제어점 그룹을 해석한다. 점 4개는 하나의 3차 Bezier 곡선이고,
    % 3n+1개 점은 이어진 3차 곡선들이다. 그 외 형태는 샘플링된 폴리라인
    % 으로 대체한다.
    curves = [];
    n = size(controlXY, 1);

    if n == 4
        curves = makeBezierCurve(cfg, id, controlXY);
    elseif mod(n - 1, 3) == 0
        for k = 1:3:(n - 3)
            curves = appendCurves(curves, makeBezierCurve(cfg, id, controlXY(k:k+3, :)));
        end
    else
        curves = makePolylineCurve(cfg, id, controlXY);
    end
end

function curve = makeBezierCurve(cfg, id, controlXY)
    % 탐색이 네 개의 제어점만이 아니라 전체 곡선 형상을 비교하도록
    % 3차 Bezier 곡선을 촘촘히 샘플링한다.
    nSamples = chooseSampleCount(cfg, controlXY);
    t = linspace(0, 1, nSamples).';
    omt = 1 - t;

    sampleXY = ...
        (omt.^3) .* controlXY(1, :) + ...
        (3 * omt.^2 .* t) .* controlXY(2, :) + ...
        (3 * omt .* t.^2) .* controlXY(3, :) + ...
        (t.^3) .* controlXY(4, :);

    curve = struct( ...
        'id', id, ...
        'type', 'cubicBezier', ...
        'controlXY', controlXY, ...
        'sampleXY', sampleXY);
end

function curve = makePolylineCurve(cfg, id, controlXY)
    % Bezier로 해석하기 어려운 제어점 배치를 위한 대체 처리다. Bezier와
    % 비슷한 간격으로 구간별 선형 경로를 샘플링한다.
    spacing = getScalar(cfg, 'alignCurveSampleSpacing', getScalar(cfg, 'sampleSpacing', 0.5));
    segLen = sqrt(sum(diff(controlXY, 1, 1).^2, 2));
    cumLen = [0; cumsum(segLen)];
    totalLen = cumLen(end);

    if totalLen <= eps
        sampleXY = controlXY(1, :);
    else
        nSamples = clampSampleCount(cfg, ceil(totalLen / spacing) + 1);
        queryLen = linspace(0, totalLen, nSamples).';
        sampleXY = [interp1(cumLen, controlXY(:, 1), queryLen, 'linear'), ...
                    interp1(cumLen, controlXY(:, 2), queryLen, 'linear')];
    end

    curve = struct( ...
        'id', id, ...
        'type', 'polylineFallback', ...
        'controlXY', controlXY, ...
        'sampleXY', sampleXY);
end

function nSamples = chooseSampleCount(cfg, controlXY)
    % 제어 다각선 길이로 곡선 길이를 근사하고 이를 샘플 개수로 변환한다.
    spacing = getScalar(cfg, 'alignCurveSampleSpacing', getScalar(cfg, 'sampleSpacing', 0.5));
    controlLen = sum(sqrt(sum(diff(controlXY, 1, 1).^2, 2)));
    nSamples = clampSampleCount(cfg, ceil(controlLen / spacing) + 1);
end

function nSamples = clampSampleCount(cfg, nSamples)
    % 너무 짧거나 긴 곡선도 수치적으로 유용하게 유지하면서 격자 탐색 비용이
    % 과도하게 커지지 않도록 샘플 개수를 제한한다.
    minSamples = getScalar(cfg, 'alignMinCurveSamples', 20);
    maxSamples = getScalar(cfg, 'alignMaxCurveSamples', 100);
    nSamples = max(minSamples, min(maxSamples, nSamples));
end

function best = gridSearch(sourceCurves, mapXY, xVals, yVals, yawVals, d2Clip, trimFraction)
    % 요각, x, y 후보를 완전 탐색하고 클리핑/절사된 곡선 집합 비용이
    % 가장 낮은 자세를 보관한다.
    best.x = [0, 0, 0];
    best.cost = inf;
    best.curveCosts = [];

    for yaw = yawVals
        rotatedCurves = rotateCurves(sourceCurves, yaw);

        for ix = 1:numel(xVals)
            x = xVals(ix);
            for iy = 1:numel(yVals)
                y = yVals(iy);
                [cost, curveCosts] = curveSetCost(rotatedCurves, mapXY, [x, y], d2Clip, trimFraction);

                if cost < best.cost
                    best.cost = cost;
                    best.x = [x, y, yaw];
                    best.curveCosts = curveCosts;
                end
            end
        end
    end
end

function [cost, curveCosts] = curveSetCost(curves, mapXY, translation, d2Clip, trimFraction)
    % 변환된 각 원본 곡선에 대해, 지도 최근접 점까지의 클리핑된 제곱거리
    % 평균을 비용으로 계산한다.
    curveCosts = zeros(numel(curves), 1);

    for i = 1:numel(curves)
        queryXY = curves(i).sampleXY + translation;
        d2 = nearestSquaredDistances(queryXY, mapXY);
        curveCosts(i) = mean(min(d2, d2Clip));
    end

    sortedCosts = sort(curveCosts);
    % 강건 집계: 일부 측정 차선이 잘린 지도에 없을 수 있으므로
    % 비용이 낮은 곡선 일부만 유지한다.
    keepCount = max(1, floor(numel(sortedCosts) * max(0, min(1, trimFraction))));
    cost = mean(sortedCosts(1:keepCount));
end

function rotatedCurves = rotateCurves(curves, yaw)
    % 요각 후보마다 곡선을 한 번 미리 회전해두면, 내부 격자에서는 병진 이동만
    % 더하면 된다.
    R = rotationMatrix(yaw);
    rotatedCurves = curves;

    for i = 1:numel(curves)
        rotatedCurves(i).sampleXY = curves(i).sampleXY * R';
        rotatedCurves(i).controlXY = curves(i).controlXY * R';
    end
end

function d2 = nearestSquaredDistances(queryXY, mapXY)
    % 질의 점에서 지도 점군까지의 최근접 이웃 제곱거리를 완전 탐색 방식으로
    % 계산한다.
    d2 = inf(size(queryXY, 1), 1);
    mapX = mapXY(:, 1);
    mapY = mapXY(:, 2);

    for i = 1:size(queryXY, 1)
        dx = mapX - queryXY(i, 1);
        dy = mapY - queryXY(i, 2);
        d2(i) = min(dx .* dx + dy .* dy);
    end
end

function outXY = transformPoints(sourceXY, x)
    % 차체 좌표계 점에 SE(2) 변환 [x, y, yaw]를 적용한다.
    R = rotationMatrix(x(3));
    outXY = sourceXY * R' + x(1:2).';
end

function R = rotationMatrix(yaw)
    % 행 벡터 형태의 점 변환에 사용할 2-D 회전 행렬이다.
    R = [cos(yaw), -sin(yaw);
         sin(yaw),  cos(yaw)];
end

function mapCrop = cropMap(mapXY, sourceXY, prior, searchX, searchY, margin)
    % 전체 탐색 영역을 포함하는 보수적인 축 정렬 자르기 영역이다. 원본 반경은
    % 탐색 중 가능한 요각 변화를 고려한다.
    sourceRadius = max(sqrt(sum(sourceXY.^2, 2)));
    xMin = prior(1) - searchX - sourceRadius - margin;
    xMax = prior(1) + searchX + sourceRadius + margin;
    yMin = prior(2) - searchY - sourceRadius - margin;
    yMax = prior(2) + searchY + sourceRadius + margin;

    keep = mapXY(:, 1) >= xMin & mapXY(:, 1) <= xMax & ...
           mapXY(:, 2) >= yMin & mapXY(:, 2) <= yMax;
    mapCrop = mapXY(keep, :);
end

function xy = collectCurveSamples(curves)
    % 최종 진단을 위해 모든 곡선 샘플점을 하나의 배열로 펼친다.
    xy = zeros(0, 2);
    for i = 1:numel(curves)
        xy = [xy; curves(i).sampleXY]; %#ok<AGROW>
    end
end

function xy = collectCurveControls(curves)
    % 그리기 함수가 각 곡선을 따로 그릴 수 있도록 NaN 구분자를 넣어
    % 제어점을 하나의 배열로 펼친다.
    xy = zeros(0, 2);
    for i = 1:numel(curves)
        xy = [xy; curves(i).controlXY; NaN, NaN]; %#ok<AGROW>
    end
end

function curves = appendCurves(curves, newCurves)
    % 파서 전반에서 사용하는 빈 배열 관례를 유지하면서 구조체 배열을
    % 이어 붙인다.
    if isempty(newCurves)
        return;
    end

    if isempty(curves)
        curves = newCurves;
    else
        curves = [curves, newCurves]; %#ok<AGROW>
    end
end

function values = makeRange(center, radius, step)
    % 중심값 주변의 양 끝을 포함하는 격자를 만든다.
    values = center + (-radius:step:radius);
    if isempty(values)
        values = center;
    end
end

function value = getScalar(s, name, defaultValue)
    % 선택적 스칼라 형태 설정값을 읽고, 없으면 기본값을 사용한다.
    value = defaultValue;
    if isfield(s, name) && ~isempty(s.(name))
        value = s.(name);
    end
end

function angle = wrapAngle(angle)
    % 요각을 [-pi, pi] 범위로 정규화한다.
    angle = atan2(sin(angle), cos(angle));
end
