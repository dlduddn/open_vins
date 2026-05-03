function [x0, info] = estimateInitialAlignment(cfg, mapDB, meas)
%ESTIMATEINITIALPOSEBYCURVELIKELIHOOD Initial SE(2) from Bezier-map likelihood.
%
% This initializer searches pose candidates and chooses the state with the
% highest curve-map likelihood. It shares the exact same likelihood model
% used later by the particle filter.

    prior = [0; 0; 0];
    x0 = prior;
    info = struct('frameIdx', [], 'x0', x0, 'prior', prior);
    
    % =========================== Integrity Check =========================== 
    [frameIdx, frame] = firstUsableBezierFrame(cfg, meas);
    if isempty(frameIdx)
        warning('estimateInitialPoseByCurveLikelihood:NoUsableQuery', ...
            'No query Bezier curve survived data association. Returning prior.');
        return;
    end

    sourceXY = collectCurveSamples(frame);
    sourceControlXY = collectCurveControls(frame);
    if isempty(sourceXY)
        warning('estimateInitialPoseByCurveLikelihood:EmptySource', ...
            'Usable query frame has no sampled curve points. Returning prior.');
        return;
    end
    % ======================================================================= 
    
    % ============================= Parameters ==============================
    searchX = getScalar(cfg, 'initSearchX', 10.0);
    searchY = getScalar(cfg, 'initSearchY', 10.0);
    searchYaw = deg2rad(getScalar(cfg, 'initSearchYawDeg', 30.0));

    fineRadiusXY = getScalar(cfg, 'initFineSerachXY', 1.0);
    fineRadiusYaw = deg2rad(getScalar(cfg, 'initFineSearchYawDeg', 2.0));
    
    cropMargin = getScalar(cfg, 'initMapCropMargin', 20.0);
    mapCrop = cropMapForSearch(mapDB, sourceXY, prior, searchX, searchY, cropMargin);
    if isempty(mapCrop), mapCrop = mapDB; end
    mapIndex = prepareCurveMapIndex(cfg, mapCrop);
    % ======================================================================= 

    % ============================= Alignment ==============================
    % Coarse
    coarseX = makeRange(prior(1), searchX, getScalar(cfg, 'initCoarseStepXY', 1.0));
    coarseY = makeRange(prior(2), searchY, getScalar(cfg, 'initCoarseStepXY', 1.0));
    coarseYaw = makeRange(prior(3), searchYaw, deg2rad(getScalar(cfg, 'initCoarseStepYawDeg', 2.0)));
    coarseBest = searchGrid(cfg, mapIndex, frame, coarseX, coarseY, coarseYaw);
    
    % Fine
    fineX = makeRange(coarseBest.x(1), fineRadiusXY, getScalar(cfg, 'initFineStepXY', 0.25));
    fineY = makeRange(coarseBest.x(2), fineRadiusXY, getScalar(cfg, 'initFineStepXY', 0.25));
    fineYaw = makeRange(coarseBest.x(3), fineRadiusYaw, deg2rad(getScalar(cfg, 'initFineStepYawDeg', 0.5)));
    if getLogical(cfg, 'initUseFineSearch', true)
        best = searchGrid(cfg, mapIndex, frame, fineX, fineY, fineYaw);
    else
        best = coarseBest;
    end

        % Results
    if ~isfinite(best.logL)
        warning('estimateInitialPoseByCurveLikelihood:NoFiniteLikelihood', ...
            'All initialization candidates failed map association. Returning prior.');
        x0 = prior;
    else
        x0 = [best.x(1); best.x(2); wrapAngle(best.x(3))];
    end

    [~, bestEval] = curveMapLikelihood(cfg, mapIndex, frame, x0);
    % ======================================================================== 

    % ================================ Info ================================ 
    info.frameIdx = frameIdx;
    info.x0 = x0;
    info.prior = prior;
    info.logLikelihood = best.logL;
    info.coarseX0 = coarseBest.x(:);
    info.coarseLogLikelihood = coarseBest.logL;
    info.numSourceCurves = numel(frame.usableIdx);
    info.numSourcePoints = size(sourceXY, 1);
    info.numMapPoints = size(mapCrop, 1);
    info.sourceCurves = frame.curves(frame.usableIdx);
    info.sourceXY = sourceXY;
    info.sourceControlXY = sourceControlXY;
    info.mapCrop = mapCrop(:, 1:2);
    info.priorAlignedXY = transformPoints(sourceXY, prior);
    info.coarseAlignedXY = transformPoints(sourceXY, coarseBest.x(:));
    info.alignedXY = transformPoints(sourceXY, x0);
    info.likelihoodInfo = bestEval;

    if isfield(bestEval, 'meanResidual') && ~isempty(bestEval.meanResidual)
        info.meanResidual = bestEval.meanResidual(1);
    else
        info.meanResidual = NaN;
    end

    fprintf(['Initial Bezier likelihood: frame=%d, x0=[%.3f %.3f %.3fdeg], ' ...
             'logL=%.2f, matchedCurves=%d/%d, samples=%d\n'], ...
        frameIdx, x0(1), x0(2), rad2deg(x0(3)), best.logL, ...
        bestEval.matchedCurveCount(1), numel(frame.usableIdx), size(sourceXY, 1));
    % ======================================================================== 
end

function [frameIdx, frame] = firstUsableBezierFrame(cfg, meas)
    frameIdx = [];
    frame = [];

    for k = 1:numel(meas)
        frameCandidate = buildQueryBezierFrame(cfg, meas(k));
        if frameCandidate.hasUsableCurves
            frameIdx = k;
            frame = frameCandidate;
            return;
        end
    end
end

function best = searchGrid(cfg, mapDB, frame, xVals, yVals, yawVals)
    best.x = [0, 0, 0];
    best.logL = -inf;

    [gridX, gridY] = ndgrid(xVals, yVals);
    xyCount = numel(gridX);
    batchSize = 512;

    for iyaw = 1:numel(yawVals)
        yaw = yawVals(iyaw);
        states = [gridX(:).'; gridY(:).'; repmat(yaw, 1, xyCount)];

        for first = 1:batchSize:xyCount
            last = min(first + batchSize - 1, xyCount);
            batchStates = states(:, first:last);
            logL = curveMapLikelihood(cfg, mapDB, frame, batchStates);
            [candidateLogL, idx] = max(logL);
            if candidateLogL > best.logL
                best.logL = candidateLogL;
                best.x = batchStates(:, idx).';
            end
        end
    end
end

function mapCrop = cropMapForSearch(mapDB, sourceXY, prior, searchX, searchY, margin)
    sourceRadius = max(sqrt(sum(sourceXY.^2, 2)));
    xMin = prior(1) - searchX - sourceRadius - margin;
    xMax = prior(1) + searchX + sourceRadius + margin;
    yMin = prior(2) - searchY - sourceRadius - margin;
    yMax = prior(2) + searchY + sourceRadius + margin;

    keep = mapDB(:, 1) >= xMin & mapDB(:, 1) <= xMax & ...
           mapDB(:, 2) >= yMin & mapDB(:, 2) <= yMax & ...
           all(isfinite(mapDB(:, 1:2)), 2);
    mapCrop = mapDB(keep, :);
end

function xy = collectCurveSamples(frame)
    xy = zeros(0, 2);
    for i = frame.usableIdx(:).'
        xy = [xy; frame.curves(i).sampleXY]; %#ok<AGROW>
    end
end

function xy = collectCurveControls(frame)
    xy = zeros(0, 2);
    for i = frame.usableIdx(:).'
        xy = [xy; frame.curves(i).controlXY; NaN, NaN]; %#ok<AGROW>
    end
end

function outXY = transformPoints(sourceXY, x)
    c = cos(x(3));
    s = sin(x(3));
    outXY = sourceXY * [c, s; -s, c] + x(1:2).';
end

function values = makeRange(center, radius, step)
    values = center + (-radius:step:radius);
    if isempty(values)
        values = center;
    end
end
