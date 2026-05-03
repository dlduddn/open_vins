function [xhat, N_eff, pfDiag] = sirLikelihood(cfg, mapDB, x0, P0, dSE2, meas, sqrtQ, len, N)
%SIRLIKELIHOOD Toy SIR particle filter using ordinary likelihood.
%
% This file intentionally does not call curveMapLikelihood().
% Measurement weights are updated directly as:
%
%     w_new = w_old .* likelihood
%
% where likelihood is computed by an analytic Bezier-to-local-line model
% implemented in this file.

    rng(cfg.pfRandomSeed);

    stateDim = 3;
    xhat = zeros(stateDim, len);
    N_eff = zeros(1, len);
    pfDiag = initPfDiagnostics(len);

    mapIndex = prepareCurveMapIndex(cfg, mapDB);

    xi = repmat(x0(:), 1, N) + sqrt(P0) * randn(stateDim, N);
    xi(3, :) = wrapAngle(xi(3, :));
    wi = (1 / N) * ones(1, N);

    for t = 1:len
        
        if t > 1
            xi = propagateParticles(xi, dSE2(:, t - 1));
            sqrtQStep = processNoiseSqrt(cfg, sqrtQ, dSE2(:, t - 1), stateDim);
            xi = addBodyFrameNoise(cfg, xi, sqrtQStep, N);
        end

        [wi, used, usableCount, meanLikelihood, matchViz] = updateWeightsLikelihood(cfg, mapIndex, meas, t, xi, wi);

        [xhat(:, t), N_eff(t)] = weightedEstimate(xi, wi);

        resampled = N_eff(t) < getScalar(cfg, 'pfResampleRatio', 0.5) * N;
        if resampled
            idx = sysresample(wi);
            xi = xi(:, idx);
            xi = roughenParticles(cfg, xi, N);
            wi = (1 / N) * ones(1, N);
        end

        visualizePfMapMatching(cfg, mapIndex, t, len, xi, wi, xhat(:, t), ...
            N_eff(t), used, resampled, usableCount, meanLikelihood, matchViz);
        pfDiag = recordPfDiagnostics(pfDiag, t, used, resampled, usableCount, ...
            meanLikelihood, N_eff(t), N, matchViz);
        logFrameProgress(cfg, t, len, used, resampled, usableCount, ...
            meanLikelihood, N_eff(t), N, xhat(:, t));
    end
end

function [wi, used, usableCount, meanLikelihood, matchViz] = ...
        updateWeightsLikelihood(cfg, mapDB, meas, t, xi, wi)
    used = false;
    usableCount = 0;
    meanLikelihood = NaN;
    matchViz = emptyMatchViz();

    [hasFrame, measFrame] = getMeasurementFrame(meas, t);
    if ~hasFrame
        return;
    end

    frame = buildQueryBezierFrame(cfg, measFrame);
    matchViz.frame = frame;
    usableCount = numel(frame.usableIdx);
    if ~frame.hasUsableCurves
        return;
    end

    [likelihood, info] = curveMapLikelihoodAnalytic(cfg, mapDB, frame, xi);
    likelihood(~isfinite(likelihood) | likelihood < 0) = 0;

    matchViz.likelihood = likelihood;
    matchViz.safeLikelihood = likelihood;
    matchViz.likelihoodInfo = info;
    if isfield(info, 'logLikelihood')
        matchViz.safeLogL = info.logLikelihood;
    end

    associatedRatio = particleAssociationRatio(cfg, info);
    matchViz.associatedRatio = associatedRatio;
    minAssociatedRatio = getScalar(cfg, 'pfMinAssociatedParticleRatio', 0.0);
    useAssociationGate = isfield(info, 'useAssociationGate') && info.useAssociationGate;
    if useAssociationGate && isfinite(associatedRatio) && associatedRatio < minAssociatedRatio
        return;
    end

    if all(likelihood <= 0)
        return;
    end

    meanLikelihood = mean(likelihood);

    wi = wi .* likelihood;
    wsum = sum(wi);
    if wsum <= 0 || ~isfinite(wsum)
        wi = (1 / numel(wi)) * ones(size(wi));
        return;
    end

    wi = wi / wsum;
    used = true;
end

function [hasFrame, measFrame] = getMeasurementFrame(meas, t)
    hasFrame = false;
    measFrame = [];

    if isstruct(meas) && isfield(meas, 'Bezier')
        if t > numel(meas.Bezier)
            return;
        end

        measFrame = meas.Bezier(t);
        if isfield(meas, 'Timestamp') && isfield(measFrame, 't_cam')
            if abs(meas.Timestamp(t) - measFrame.t_cam) > eps(max(abs(meas.Timestamp(t)), 1))
                measFrame = [];
                return;
            end
        end

        hasFrame = true;
        return;
    end

    if t <= numel(meas)
        measFrame = meas(t);
        hasFrame = true;
    end
end

function [likelihood, info] = curveMapLikelihoodAnalytic(cfg, mapDB, queryInput, xStates)
    if isvector(xStates)
        xStates = xStates(:);
    end
    if size(xStates, 1) ~= 3 && size(xStates, 2) == 3
        xStates = xStates.';
    end

    nState = size(xStates, 2);
    likelihood = zeros(1, nState);
    info = emptyAnalyticInfo(nState);

    mapIndex = asAnalyticMapIndex(cfg, mapDB);
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
    minMatchedCurves = getScalar(cfg, 'assocMinMatchedCurves', 1);
    maxEffSamples = getScalar(cfg, 'likelihoodMaxEffectiveSamples', inf);
    missPenalty = -abs(getScalar(cfg, 'likelihoodMissLogPenalty', 60.0));
    segmentFitK = max(2, round(getScalar(cfg, 'assocMapSegmentFitK', 5)));
    useAssociationGate = getLogical(cfg, 'assocUseCurveMapCorrespondenceRule', false);
    useFovGate = getLogical(cfg, 'assocUseCameraFovGate', false);
    if useFovGate
        fovGate = computeCameraFovGate(cfg);
    else
        fovGate = disabledAnalyticFovGate();
    end
    info.useAssociationGate = useAssociationGate;
    info.useFovGate = fovGate.enabled;

    finiteState = all(isfinite(xStates), 1);
    evalIdx = find(finiteState);
    if isempty(evalIdx)
        return;
    end

    xEval = xStates(:, evalIdx);
    c = cos(xEval(3, :));
    s = sin(xEval(3, :));

    logScore = zeros(1, nState);
    matchedCount = zeros(1, nState);
    residualSum = zeros(1, nState);
    residualSquaredSum = zeros(1, nState);
    residualCount = zeros(1, nState);

    for i = 1:numel(curves)
        controlXY = getCurveControlXY(curves(i));
        if isempty(controlXY)
            continue;
        end
        if fovGate.enabled && ~analyticCurveVisibleInBody(controlXY, fovGate)
            continue;
        end

        coeff = bezierPolynomialCoefficients(controlXY);
        anchorBody = bezierPointAt(controlXY, 0.5);
        anchorMap = transformBodyPointBatch(anchorBody, xEval, c, s);

        [lineCenter, lineNormal, nnIdx] = fitLocalMapLines(mapIndex, anchorMap, segmentFitK);
        validLine = all(isfinite(lineCenter), 2).' & all(isfinite(lineNormal), 2).';
        if fovGate.enabled
            fovKeep = nearestAnalyticLinesInsideFov(mapIndex, xEval, nnIdx, ...
                lineCenter, fovGate, maxDist);
            validLine = validLine & fovKeep(:).';
        end
        validLine = validLine(:).';
        if ~any(validLine)
            continue;
        end

        meanSquaredResidual = analyticBezierLineMeanSquare(coeff, xEval, c, s, ...
            lineCenter, lineNormal);
        meanSquaredResidual = meanSquaredResidual(:).';
        rmsResidual = sqrt(meanSquaredResidual);

        if useAssociationGate
            scored = validLine & rmsResidual <= maxDist;
        else
            scored = validLine;
        end
        scored = scored(:).';
        if ~any(scored)
            continue;
        end

        usedEvalIdx = evalIdx(scored);
        nEff = analyticCurveEffectiveSampleCount(cfg, curves(i), maxEffSamples);
        curveLogScore = -0.5 * nEff .* meanSquaredResidual(scored) ./ (sigma ^ 2);

        logScore(usedEvalIdx) = logScore(usedEvalIdx) + curveLogScore;
        matchedCount(usedEvalIdx) = matchedCount(usedEvalIdx) + 1;
        residualSquaredSum(usedEvalIdx) = residualSquaredSum(usedEvalIdx) + meanSquaredResidual(scored);
        residualSum(usedEvalIdx) = residualSum(usedEvalIdx) + rmsResidual(scored);
        residualCount(usedEvalIdx) = residualCount(usedEvalIdx) + 1;
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
    if ~any(valid)
        logScore(~finiteState) = -inf;
        info.logLikelihood = logScore;
        return;
    end

    rejected = finiteState & ~valid;
    logScore(rejected) = missPenalty;
    logScore(~finiteState) = -inf;
    info.logLikelihood = logScore;
    likelihood = logScoreToRelativeLikelihood(logScore);
end

function mapIndex = asAnalyticMapIndex(cfg, mapDB)
    if isstruct(mapDB) && isfield(mapDB, 'xy')
        mapIndex = mapDB;
    else
        mapIndex = prepareCurveMapIndex(cfg, mapDB);
    end

    if ~isfield(mapIndex, 'xy') || size(mapIndex.xy, 2) < 2
        mapIndex = struct('xy', zeros(0, 2), 'hasKdTree', false, 'searcher', []);
    end
    if ~isfield(mapIndex, 'hasKdTree')
        mapIndex.hasKdTree = false;
    end
    if ~isfield(mapIndex, 'searcher')
        mapIndex.searcher = [];
    end
end

function controlXY = getCurveControlXY(curve)
    controlXY = zeros(0, 2);
    if isfield(curve, 'controlXY') && size(curve.controlXY, 1) == 4 && size(curve.controlXY, 2) >= 2
        controlXY = curve.controlXY(:, 1:2);
    end
    if any(~isfinite(controlXY(:)))
        controlXY = zeros(0, 2);
    end
end

function coeff = bezierPolynomialCoefficients(controlXY)
    p0 = controlXY(1, :);
    p1 = controlXY(2, :);
    p2 = controlXY(3, :);
    p3 = controlXY(4, :);
    coeff = [
        p0;
        -3 * p0 + 3 * p1;
        3 * p0 - 6 * p1 + 3 * p2;
        -p0 + 3 * p1 - 3 * p2 + p3
    ];
end

function pointXY = bezierPointAt(controlXY, t)
    omt = 1 - t;
    pointXY = ...
        (omt ^ 3) * controlXY(1, :) + ...
        (3 * omt ^ 2 * t) * controlXY(2, :) + ...
        (3 * omt * t ^ 2) * controlXY(3, :) + ...
        (t ^ 3) * controlXY(4, :);
end

function mapXY = transformBodyPointBatch(bodyXY, xStates, c, s)
    mapXY = [
        bodyXY(1) .* c - bodyXY(2) .* s + xStates(1, :);
        bodyXY(1) .* s + bodyXY(2) .* c + xStates(2, :)
    ].';
end

function [lineCenter, lineNormal, nnIdx] = fitLocalMapLines(mapIndex, anchorXY, segmentFitK)
    k = min(max(2, segmentFitK), size(mapIndex.xy, 1));
    nnIdx = nearestKMapIndices(mapIndex, anchorXY, k);
    lineCenter = NaN(size(anchorXY));
    lineNormal = NaN(size(anchorXY));
    if isempty(nnIdx)
        return;
    end

    neighborX = reshape(mapIndex.xy(nnIdx(:), 1), size(nnIdx));
    neighborY = reshape(mapIndex.xy(nnIdx(:), 2), size(nnIdx));

    centerX = mean(neighborX, 2);
    centerY = mean(neighborY, 2);
    dx = neighborX - centerX;
    dy = neighborY - centerY;

    covXX = sum(dx .* dx, 2);
    covXY = sum(dx .* dy, 2);
    covYY = sum(dy .* dy, 2);
    theta = 0.5 * atan2(2.0 * covXY, covXX - covYY);

    dirX = cos(theta);
    dirY = sin(theta);
    degenerate = (covXX + covYY) <= eps;

    lineCenter = [centerX, centerY];
    lineNormal = [-dirY, dirX];
    lineNormal(degenerate, :) = NaN;
end

function nnIdx = nearestKMapIndices(mapIndex, queryXY, k)
    if isempty(queryXY) || isempty(mapIndex.xy)
        nnIdx = zeros(0, k);
        return;
    end

    if mapIndex.hasKdTree && ~isempty(mapIndex.searcher) && exist('knnsearch', 'file') == 2
        nnIdx = knnsearch(mapIndex.searcher, queryXY, 'K', k);
        nnIdx = reshape(nnIdx, size(queryXY, 1), k);
        return;
    end

    nnIdx = ones(size(queryXY, 1), k);
    mapX = mapIndex.xy(:, 1).';
    mapY = mapIndex.xy(:, 2).';
    chunkSize = 256;
    for first = 1:chunkSize:size(queryXY, 1)
        last = min(first + chunkSize - 1, size(queryXY, 1));
        dx = queryXY(first:last, 1) - mapX;
        dy = queryXY(first:last, 2) - mapY;
        d2 = dx .* dx + dy .* dy;
        [~, order] = sort(d2, 2, 'ascend');
        nnIdx(first:last, :) = order(:, 1:k);
    end
end

function meanSquare = analyticBezierLineMeanSquare(coeff, xStates, c, s, lineCenter, lineNormal)
    nState = size(xStates, 2);
    residualCoeff = zeros(4, nState);

    for k = 1:4
        body = coeff(k, :);
        mapCoeff = [
            body(1) .* c - body(2) .* s;
            body(1) .* s + body(2) .* c
        ].';
        if k == 1
            mapCoeff = mapCoeff + xStates(1:2, :).' - lineCenter;
        end
        residualCoeff(k, :) = sum(mapCoeff .* lineNormal, 2).';
    end

    meanSquare = integrateCubicSquared(residualCoeff);
end

function integralValue = integrateCubicSquared(coeff)
    integralValue = zeros(1, size(coeff, 2));
    for i = 1:4
        for j = 1:4
            integralValue = integralValue + coeff(i, :) .* coeff(j, :) ./ (i + j - 1);
        end
    end
    integralValue = max(integralValue, 0.0);
end

function nEff = analyticCurveEffectiveSampleCount(cfg, curve, maxEffSamples)
    spacing = max(getScalar(cfg, 'assocCurveSampleSpacing', getScalar(cfg, 'sampleSpacing', 0.5)), eps);
    if isfield(curve, 'length') && isfinite(curve.length) && curve.length > 0
        nEff = max(1, ceil(curve.length / spacing) + 1);
    elseif isfield(curve, 'sampleXY') && ~isempty(curve.sampleXY)
        nEff = size(curve.sampleXY, 1);
    else
        nEff = 1;
    end
    nEff = min(nEff, maxEffSamples);
end

function keep = nearestAnalyticLinesInsideFov(mapIndex, xStates, nnIdx, lineCenter, fovGate, rangeMargin)
    nearestXY = mapIndex.xy(nnIdx(:, 1), :);
    lineCenterBody = transformMapPointsToBody(lineCenter, xStates);
    nearestBody = transformMapPointsToBody(nearestXY, xStates);
    keep = pointsInsideCameraFovAnalytic(lineCenterBody, fovGate, rangeMargin) | ...
        pointsInsideCameraFovAnalytic(nearestBody, fovGate, rangeMargin);
    keep = keep(:).';
end

function bodyXY = transformMapPointsToBody(mapXY, xStates)
    c = cos(xStates(3, :)).';
    s = sin(xStates(3, :)).';
    rel = mapXY - xStates(1:2, :).';
    bodyXY = [rel(:, 1) .* c + rel(:, 2) .* s, ...
             -rel(:, 1) .* s + rel(:, 2) .* c];
end

function visible = analyticCurveVisibleInBody(controlXY, fovGate)
    probeXY = [
        controlXY(1, :);
        bezierPointAt(controlXY, 0.5);
        controlXY(4, :)
    ];
    visible = any(pointsInsideCameraFovAnalytic(probeXY, fovGate, 0.0));
end

function keep = pointsInsideCameraFovAnalytic(bodyXY, gate, rangeMargin)
    rangeMargin = max(rangeMargin, 0.0);
    relAngle = wrapAngle(atan2(bodyXY(:, 2), bodyXY(:, 1)) - gate.yawOffset);
    range = sqrt(sum(bodyXY .^ 2, 2));
    forward = bodyXY(:, 1) * cos(gate.yawOffset) + bodyXY(:, 2) * sin(gate.yawOffset);

    keep = forward > 0 & ...
        range >= max(0.0, gate.minRange - rangeMargin) & ...
        range <= gate.maxRange + rangeMargin & ...
        relAngle >= -gate.rightAngle - gate.margin & ...
        relAngle <= gate.leftAngle + gate.margin;
end

function likelihood = logScoreToRelativeLikelihood(logScore)
    likelihood = zeros(size(logScore));
    finiteMask = isfinite(logScore);
    if ~any(finiteMask)
        return;
    end

    bestScore = max(logScore(finiteMask));
    likelihood(finiteMask) = exp(logScore(finiteMask) - bestScore);
    likelihood(~isfinite(likelihood)) = 0;
end

function info = emptyAnalyticInfo(nState)
    info = struct( ...
        'hasMeasurement', false, ...
        'usableCurveCount', 0, ...
        'useAssociationGate', false, ...
        'useFovGate', false, ...
        'matchedCurveCount', zeros(1, nState), ...
        'meanResidual', inf(1, nState), ...
        'meanSquaredResidual', inf(1, nState), ...
        'logLikelihood', -inf(1, nState));
end

function gate = disabledAnalyticFovGate()
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

function xi = propagateParticles(xi, dPose)
    dX = dPose(1);
    dY = dPose(2);
    dYaw = dPose(3);

    c = cos(xi(3, :));
    s = sin(xi(3, :));
    xi(1, :) = xi(1, :) + c .* dX - s .* dY;
    xi(2, :) = xi(2, :) + s .* dX + c .* dY;
    xi(3, :) = wrapAngle(xi(3, :) + dYaw);
end

function sqrtQStep = processNoiseSqrt(cfg, sqrtQ, dPose, stateDim)
    if ~isempty(sqrtQ)
        if ~isequal(size(sqrtQ), [stateDim, stateDim])
            error('sqrtQ must be a %d x %d matrix.', stateDim, stateDim);
        end
        sqrtQStep = sqrtQ;
        return;
    end

    dX = dPose(1);
    dY = dPose(2);
    dYaw = dPose(3);
    stepDist = hypot(dX, dY);

    stdForward = getScalar(cfg, 'pfProcessForwardStd', 0.0) + ...
        getScalar(cfg, 'pfProcessTransScale', 0.0) * stepDist;
    stdLateral = getScalar(cfg, 'pfProcessLateralStd', 0.0) + ...
        getScalar(cfg, 'pfProcessTransScale', 0.0) * stepDist;
    stdYaw = deg2rad(getScalar(cfg, 'pfProcessYawStdDeg', 0.0)) + ...
        getScalar(cfg, 'pfProcessYawScale', 0.0) * abs(dYaw);

    sqrtQStep = diag(max([stdForward, stdLateral, stdYaw], 0.0));
end

function xi = addBodyFrameNoise(cfg, xi, sqrtQStep, N)
    if isempty(sqrtQStep) || ~any(sqrtQStep(:))
        return;
    end

    noise = sqrtQStep * randn(3, N);
    noiseFrame = lower(string(getField(cfg, 'pfProcessNoiseFrame', 'body')));
    if noiseFrame == "map"
        xi = xi + noise;
        xi(3, :) = wrapAngle(xi(3, :));
        return;
    end

    c = cos(xi(3, :));
    s = sin(xi(3, :));
    xi(1, :) = xi(1, :) + c .* noise(1, :) - s .* noise(2, :);
    xi(2, :) = xi(2, :) + s .* noise(1, :) + c .* noise(2, :);
    xi(3, :) = wrapAngle(xi(3, :) + noise(3, :));
end

function xi = roughenParticles(cfg, xi, N)
    if ~getLogical(cfg, 'pfRoughenAfterResample', false)
        return;
    end

    stdRough = [
        getScalar(cfg, 'pfRoughenXStd', 0.0);
        getScalar(cfg, 'pfRoughenYStd', 0.0);
        deg2rad(getScalar(cfg, 'pfRoughenYawStdDeg', 0.0))
    ];
    if ~any(stdRough > 0)
        return;
    end

    xi = xi + diag(stdRough) * randn(3, N);
    xi(3, :) = wrapAngle(xi(3, :));
end

function [x, neff] = weightedEstimate(xi, wi)
    x = zeros(3, 1);
    x(1:2) = sum(xi(1:2, :) .* wi, 2);
    x(3) = atan2(sum(sin(xi(3, :)) .* wi), sum(cos(xi(3, :)) .* wi));
    neff = 1 / sum(wi .^ 2);
end

function pfDiag = initPfDiagnostics(len)
    pfDiag = struct( ...
        'usedUpdate', false(1, len), ...
        'resampled', false(1, len), ...
        'usableCurves', zeros(1, len), ...
        'meanLikelihood', NaN(1, len), ...
        'meanLogL', NaN(1, len), ...
        'nEff', NaN(1, len), ...
        'nEffRatio', NaN(1, len), ...
        'associatedParticleRatio', NaN(1, len), ...
        'matchedCurveCountMean', NaN(1, len), ...
        'meanResidual', NaN(1, len), ...
        'meanSquaredResidual', NaN(1, len));
end

function pfDiag = recordPfDiagnostics(pfDiag, t, used, resampled, usableCount, ...
        meanLikelihood, neff, nParticle, matchViz)
    pfDiag.usedUpdate(t) = logical(used);
    pfDiag.resampled(t) = logical(resampled);
    pfDiag.usableCurves(t) = usableCount;
    pfDiag.meanLikelihood(t) = meanLikelihood;
    pfDiag.meanLogL(t) = meanLikelihood;
    pfDiag.nEff(t) = neff;
    pfDiag.nEffRatio(t) = neff / max(nParticle, 1);

    if isstruct(matchViz) && isfield(matchViz, 'associatedRatio')
        pfDiag.associatedParticleRatio(t) = matchViz.associatedRatio;
    end

    if isstruct(matchViz) && isfield(matchViz, 'likelihoodInfo') && ...
            isstruct(matchViz.likelihoodInfo)
        info = matchViz.likelihoodInfo;
        if isfield(info, 'matchedCurveCount')
            pfDiag.matchedCurveCountMean(t) = mean(info.matchedCurveCount, 'omitnan');
        end
        if isfield(info, 'meanResidual')
            pfDiag.meanResidual(t) = mean(info.meanResidual, 'omitnan');
        end
        if isfield(info, 'meanSquaredResidual')
            pfDiag.meanSquaredResidual(t) = mean(info.meanSquaredResidual, 'omitnan');
        end
    end
end

function matchViz = emptyMatchViz()
    matchViz = struct( ...
        'frame', [], ...
        'likelihood', [], ...
        'safeLikelihood', [], ...
        'likelihoodInfo', [], ...
        'associatedRatio', NaN);
end

function ratio = particleAssociationRatio(cfg, likelihoodInfo)
    ratio = NaN;
    if ~isstruct(likelihoodInfo) || ~isfield(likelihoodInfo, 'matchedCurveCount')
        return;
    end

    minMatchedCurves = getScalar(cfg, 'assocMinMatchedCurves', 1);
    matched = likelihoodInfo.matchedCurveCount >= minMatchedCurves;
    if isempty(matched)
        return;
    end

    ratio = mean(matched);
end

function logFrameProgress(cfg, t, len, used, resampled, usableCount, ...
        meanLikelihood, neff, nParticle, x)
    if ~getLogical(cfg, 'pfLogProgress', false) || ...
            ~shouldLogFrame(t, len, getScalar(cfg, 'pfLogInterval', 50))
        return;
    end

    fprintf(['PF Likelihood [%5d/%5d] used=%d curves=%d N_eff=%.1f/%d ' ...
             'meanL=%.4g x=[%.3f %.3f %.2fdeg] resampled=%d\n'], ...
        t, len, used, usableCount, neff, nParticle, meanLikelihood, ...
        x(1), x(2), rad2deg(x(3)), resampled);
end

function doLog = shouldLogFrame(t, len, logInterval)
    doLog = t == 1 || t == len || mod(t, logInterval) == 0;
end
