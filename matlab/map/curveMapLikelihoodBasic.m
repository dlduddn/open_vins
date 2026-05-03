function [likelihood, info] = curveMapLikelihoodBasic(cfg, mapDB, queryInput, xStates)
%CURVEMAPLIKELIHOODBASIC Simple Gaussian likelihood with association/FOV gates.
%
% This function deliberately does not call curveMapLikelihood().
% It keeps only the basic pieces:
%
%   1. transform query curve samples by each particle pose,
%   2. apply the camera FOV gate,
%   3. apply a curve-to-map association gate,
%   4. compute ordinary likelihood from nearest-map-point residuals.
%
% The likelihood model is:
%
%   L = exp(-0.5 * mean(d^2) / sigma^2)

    if isvector(xStates)
        xStates = xStates(:);
    end
    if size(xStates, 1) ~= 3 && size(xStates, 2) == 3
        xStates = xStates.';
    end

    nState = size(xStates, 2);
    likelihood = zeros(1, nState);
    info = emptyInfo(nState);

    mapIndex = asMapIndexBasic(cfg, mapDB);
    if isempty(mapIndex.xy)
        return;
    end

    frame = queryInput;
    if ~isstruct(frame) || ~isfield(frame, 'curves')
        frame = buildQueryBezierFrame(cfg, queryInput);
    end

    info.hasMeasurement = isfield(frame, 'hasUsableCurves') && frame.hasUsableCurves;
    if ~info.hasMeasurement
        return;
    end

    if isfield(frame, 'usableIdx') && ~isempty(frame.usableIdx)
        curves = frame.curves(frame.usableIdx);
    else
        curves = frame.curves;
    end
    info.usableCurveCount = numel(curves);
    if isempty(curves)
        return;
    end

    sigma = max(getScalar(cfg, 'likelihoodSigma', 1.0), eps);
    maxDist = getScalar(cfg, 'assocMaxCurveMapDist', 2.5);
    maxD2 = maxDist ^ 2;
    minMatchFraction = getScalar(cfg, 'assocMinCurveMatchFraction', 0.45);
    minMatchedCurves = getScalar(cfg, 'assocMinMatchedCurves', 1);
    useAssociationGate = getLogical(cfg, 'assocUseCurveMapCorrespondenceRule', false);
    useFovGate = getLogical(cfg, 'assocUseCameraFovGate', false);
    if useFovGate
        fovGate = computeCameraFovGate(cfg);
    else
        fovGate = disabledFovGate();
    end
    info.useAssociationGate = useAssociationGate;
    info.useFovGate = fovGate.enabled;
    curveSamples = prepareCurveSamples(curves, fovGate);

    finiteState = all(isfinite(xStates), 1);
    evalIdx = find(finiteState);
    if isempty(evalIdx)
        return;
    end

    matchedCount = zeros(1, nState);
    residualSum = zeros(1, nState);
    residualSquaredSum = zeros(1, nState);
    residualCount = zeros(1, nState);

    for i = 1:numel(curves)
        sampleXY = curveSamples{i};
        if isempty(sampleXY)
            continue;
        end

        nSample = size(sampleXY, 1);
        xEval = xStates(:, evalIdx);
        c = cos(xEval(3, :));
        s = sin(xEval(3, :));

        qx = sampleXY(:, 1) .* c - sampleXY(:, 2) .* s + xEval(1, :);
        qy = sampleXY(:, 1) .* s + sampleXY(:, 2) .* c + xEval(2, :);

        [d2, nnIdx] = nearestMapSquaredDistancesBasic(mapIndex, [qx(:), qy(:)]);

        if fovGate.enabled
            keepNearest = nearestPointsInsideFovBasic(mapIndex, xEval, nnIdx, nSample, fovGate, maxDist);
            d2(~keepNearest) = inf;
        end

        d2 = reshape(d2, nSample, numel(evalIdx));
        dist = sqrt(d2);
        matchFraction = mean(dist <= maxDist, 1);
        if useAssociationGate
            scored = matchFraction >= minMatchFraction;
        else
            scored = true(1, numel(evalIdx));
        end

        finiteResidual = isfinite(d2);
        curveResidualCount = sum(finiteResidual, 1);
        curveUsed = scored & curveResidualCount > 0;
        if ~any(curveUsed)
            continue;
        end

        usedEvalIdx = evalIdx(curveUsed);
        d2Used = d2(:, curveUsed);
        d2Used(~isfinite(d2Used)) = 0;
        d2Clipped = min(d2Used, maxD2);

        matchedCount(usedEvalIdx) = matchedCount(usedEvalIdx) + 1;
        residualSquaredSum(usedEvalIdx) = residualSquaredSum(usedEvalIdx) + sum(d2Clipped, 1);
        residualSum(usedEvalIdx) = residualSum(usedEvalIdx) + sum(sqrt(d2Clipped), 1);
        residualCount(usedEvalIdx) = residualCount(usedEvalIdx) + curveResidualCount(curveUsed);
    end

    info.matchedCurveCount = matchedCount;
    hasResidual = residualCount > 0;
    info.meanSquaredResidual(hasResidual) = residualSquaredSum(hasResidual) ./ residualCount(hasResidual);
    info.meanResidual(hasResidual) = residualSum(hasResidual) ./ residualCount(hasResidual);

    associated = true(1, nState);
    if useAssociationGate
        associated = matchedCount >= minMatchedCurves;
    end

    valid = hasResidual & associated;
    likelihood(valid) = exp(-0.5 * info.meanSquaredResidual(valid) / (sigma ^ 2));
end

function mapIndex = asMapIndexBasic(cfg, mapDB)
    if isstruct(mapDB) && isfield(mapDB, 'xy')
        mapIndex = mapDB;
    else
        mapIndex = prepareCurveMapIndex(cfg, mapDB);
    end

    if ~isfield(mapIndex, 'xy') || size(mapIndex.xy, 2) < 2
        mapIndex = emptyMapIndexBasic();
    end
    if ~isfield(mapIndex, 'hasKdTree')
        mapIndex.hasKdTree = false;
    end
    if ~isfield(mapIndex, 'searcher')
        mapIndex.searcher = [];
    end
end

function mapIndex = emptyMapIndexBasic()
    mapIndex = struct( ...
        'xy', zeros(0, 2), ...
        'hasKdTree', false, ...
        'searcher', []);
end

function curveSamples = prepareCurveSamples(curves, fovGate)
    curveSamples = cell(1, numel(curves));
    for i = 1:numel(curves)
        sampleXY = getCurveSamples(curves(i));
        if fovGate.enabled && ~isempty(sampleXY)
            sampleXY = sampleXY(pointsInsideCameraFov(sampleXY, fovGate, 0.0), :);
        end
        curveSamples{i} = sampleXY;
    end
end

function sampleXY = getCurveSamples(curve)
    sampleXY = zeros(0, 2);
    if ~isfield(curve, 'sampleXY') || size(curve.sampleXY, 2) < 2
        return;
    end

    sampleXY = curve.sampleXY(:, 1:2);
    sampleXY = sampleXY(all(isfinite(sampleXY), 2), :);
end

function [d2, idx] = nearestMapSquaredDistancesBasic(mapIndex, queryXY)
    d2 = inf(size(queryXY, 1), 1);
    idx = ones(size(queryXY, 1), 1);
    if isempty(mapIndex.xy) || isempty(queryXY)
        return;
    end

    if mapIndex.hasKdTree && ~isempty(mapIndex.searcher) && exist('knnsearch', 'file') == 2
        [idx, dist] = knnsearch(mapIndex.searcher, queryXY, 'K', 1);
        d2 = dist .^ 2;
        return;
    end

    mapX = mapIndex.xy(:, 1).';
    mapY = mapIndex.xy(:, 2).';
    chunkSize = 256;

    for first = 1:chunkSize:size(queryXY, 1)
        last = min(first + chunkSize - 1, size(queryXY, 1));
        dx = queryXY(first:last, 1) - mapX;
        dy = queryXY(first:last, 2) - mapY;
        [d2(first:last), localIdx] = min(dx .* dx + dy .* dy, [], 2);
        idx(first:last) = localIdx;
    end
end

function keep = nearestPointsInsideFovBasic(mapIndex, xStates, nnIdx, nSample, fovGate, rangeMargin)
    nState = size(xStates, 2);
    stateIdx = repelem(1:nState, nSample).';
    nearestXY = mapIndex.xy(nnIdx, :);
    stateXY = xStates(1:2, stateIdx).';
    yaw = xStates(3, stateIdx).';

    rel = nearestXY - stateXY;
    c = cos(yaw);
    s = sin(yaw);
    bodyXY = [rel(:, 1) .* c + rel(:, 2) .* s, ...
             -rel(:, 1) .* s + rel(:, 2) .* c];
    keep = pointsInsideCameraFov(bodyXY, fovGate, rangeMargin);
end

function keep = pointsInsideCameraFov(bodyXY, gate, rangeMargin)
    rangeMargin = max(rangeMargin, 0.0);
    relAngle = wrapAngle(atan2(bodyXY(:, 2), bodyXY(:, 1)) - gate.yawOffset);
    range = sqrt(sum(bodyXY .^ 2, 2));
    forward = bodyXY(:, 1) * cos(gate.yawOffset) + ...
        bodyXY(:, 2) * sin(gate.yawOffset);

    keep = forward > 0 & ...
        range >= max(0.0, gate.minRange - rangeMargin) & ...
        range <= gate.maxRange + rangeMargin & ...
        relAngle >= -gate.rightAngle - gate.margin & ...
        relAngle <= gate.leftAngle + gate.margin;
end

function info = emptyInfo(nState)
    info = struct( ...
        'hasMeasurement', false, ...
        'usableCurveCount', 0, ...
        'useAssociationGate', false, ...
        'useFovGate', false, ...
        'matchedCurveCount', zeros(1, nState), ...
        'meanResidual', inf(1, nState), ...
        'meanSquaredResidual', inf(1, nState));
end

function gate = disabledFovGate()
    gate = struct( ...
        'enabled', false, ...
        'leftAngle', NaN, ...
        'rightAngle', NaN, ...
        'yawOffset', 0.0, ...
        'margin', 0.0, ...
        'minRange', 0.0, ...
        'maxRange', inf, ...
        'mode', 'nearest');
end
