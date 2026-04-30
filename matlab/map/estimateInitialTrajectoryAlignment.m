function [x0, info] = estimateInitialTrajectoryAlignment(cfg, mapDB, meas, dSE2, len)
%ESTIMATEINITIALTRAJECTORYALIGNMENT Initial SE(2) from several early frames.
%
% Single-frame curve matching is weak on straight roads because translation
% along the lane direction is poorly observable. This initializer scores an
% initial pose by propagating it through early VIO odometry and accumulating
% curve-map likelihood over several usable frames. A weak prior keeps the
% unobservable longitudinal component from overfitting inference artifacts.

    prior = [
        getScalar(cfg, 'initPriorX', 0.0);
        getScalar(cfg, 'initPriorY', 0.0);
        deg2rad(getScalar(cfg, 'initPriorYawDeg', 0.0))
    ];
    x0 = prior;

    info = struct('frameIdx', 1, 'x0', x0, 'prior', prior, ...
        'selectedFrameIdx', [], 'logLikelihood', -inf);

    if nargin < 5 || isempty(len)
        len = numel(meas);
    end

    relPose = relativeTrajectory(dSE2, len);
    [frameIdx, frames] = selectInitialFrames(cfg, meas, len);
    if isempty(frameIdx)
        warning('estimateInitialTrajectoryAlignment:NoUsableQuery', ...
            'No query Bezier curve survived data association. Returning prior.');
        return;
    end

    mapCrop = cropMapForTrajectorySearch(cfg, mapDB, frames, frameIdx, relPose, prior);
    if isempty(mapCrop)
        mapCrop = mapDB;
    end
    mapIndex = prepareCurveMapIndex(cfg, mapCrop);

    searchX = getScalar(cfg, 'initSearchX', 10.0);
    searchY = getScalar(cfg, 'initSearchY', 10.0);
    searchYaw = deg2rad(getScalar(cfg, 'initSearchYawDeg', 30.0));

    coarseX = makeRange(prior(1), searchX, getScalar(cfg, 'initCoarseStepXY', 1.0));
    coarseY = makeRange(prior(2), searchY, getScalar(cfg, 'initCoarseStepXY', 1.0));
    coarseYaw = makeRange(prior(3), searchYaw, deg2rad(getScalar(cfg, 'initCoarseStepYawDeg', 2.0)));
    coarseBest = searchInitialGrid(cfg, mapIndex, frames, frameIdx, relPose, coarseX, coarseY, coarseYaw, prior);

    fineRadiusXY = getScalar(cfg, 'initFineRadiusXY', 1.0);
    fineRadiusYaw = deg2rad(getScalar(cfg, 'initFineRadiusYawDeg', 2.0));
    fineX = makeRange(coarseBest.x(1), fineRadiusXY, getScalar(cfg, 'initFineStepXY', 0.25));
    fineY = makeRange(coarseBest.x(2), fineRadiusXY, getScalar(cfg, 'initFineStepXY', 0.25));
    fineYaw = makeRange(coarseBest.x(3), fineRadiusYaw, deg2rad(getScalar(cfg, 'initFineStepYawDeg', 0.5)));
    best = searchInitialGrid(cfg, mapIndex, frames, frameIdx, relPose, fineX, fineY, fineYaw, prior);

    if ~isfinite(best.logL)
        warning('estimateInitialTrajectoryAlignment:NoFiniteLikelihood', ...
            'All initialization candidates failed map association. Returning prior.');
        x0 = prior;
    else
        x0 = [best.x(1); best.x(2); wrapAngle(best.x(3))];
    end

    sourceXY = collectTrajectorySourceXY(frames, frameIdx, relPose);
    sourceControlXY = collectTrajectoryControlXY(frames, frameIdx, relPose);

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

    evalInfo = evaluateSelectedFrames(cfg, mapIndex, frames, frameIdx, relPose, x0);
    info.likelihoodInfo = evalInfo;
    info.meanResidual = evalInfo.meanResidual;

    fprintf(['Initial trajectory likelihood: frames=%d, idx=[%d..%d], ' ...
             'x0=[%.3f %.3f %.3fdeg], logL=%.2f, curves=%d, samples=%d\n'], ...
        numel(frameIdx), frameIdx(1), frameIdx(end), x0(1), x0(2), ...
        rad2deg(x0(3)), best.logL, info.numSourceCurves, size(sourceXY, 1));
end

function relPose = relativeTrajectory(dSE2, len)
    relPose = zeros(3, len);
    if isempty(dSE2)
        return;
    end

    for t = 2:len
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
        pick = unique(round(linspace(1, numel(usableIdx), maxFrames)));
        selectedIdx = usableIdx(pick);
        selectedFrames = frameCache(pick);
    else
        selectedIdx = usableIdx;
        selectedFrames = frameCache;
    end
end

function best = searchInitialGrid(cfg, mapIndex, frames, frameIdx, relPose, xVals, yVals, yawVals, prior)
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
    stdX = max(getScalar(cfg, 'initPriorStdX', inf), eps);
    stdY = max(getScalar(cfg, 'initPriorStdY', inf), eps);
    stdYaw = max(deg2rad(getScalar(cfg, 'initPriorStdYawDeg', inf)), eps);

    dx = states(1, :) - prior(1);
    dy = states(2, :) - prior(2);
    dyaw = wrapAngle(states(3, :) - prior(3));
    score = -0.5 * ((dx / stdX) .^ 2 + (dy / stdY) .^ 2 + (dyaw / stdYaw) .^ 2);
end

function statesFrame = propagateInitialBatch(states0, rel)
    c = cos(states0(3, :));
    s = sin(states0(3, :));

    statesFrame = states0;
    statesFrame(1, :) = states0(1, :) + c * rel(1) - s * rel(2);
    statesFrame(2, :) = states0(2, :) + s * rel(1) + c * rel(2);
    statesFrame(3, :) = wrapAngle(states0(3, :) + rel(3));
end

function evalInfo = evaluateSelectedFrames(cfg, mapIndex, frames, frameIdx, relPose, x0)
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
    sourceXY = zeros(0, 2);
    for j = 1:numel(frames)
        rel = relPose(:, frameIdx(j));
        for i = frames(j).usableIdx(:).'
            sourceXY = [sourceXY; transformBodyToInitial(frames(j).curves(i).sampleXY, rel)]; %#ok<AGROW>
        end
    end
end

function sourceXY = collectTrajectoryControlXY(frames, frameIdx, relPose)
    sourceXY = zeros(0, 2);
    for j = 1:numel(frames)
        rel = relPose(:, frameIdx(j));
        for i = frames(j).usableIdx(:).'
            sourceXY = [sourceXY; transformBodyToInitial(frames(j).curves(i).controlXY, rel); NaN, NaN]; %#ok<AGROW>
        end
    end
end

function outXY = transformBodyToInitial(bodyXY, rel)
    c = cos(rel(3));
    s = sin(rel(3));
    outXY = bodyXY * [c, s; -s, c] + rel(1:2).';
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

function value = getScalar(s, name, defaultValue)
    value = defaultValue;
    if isfield(s, name) && ~isempty(s.(name))
        value = s.(name);
    end
end

function angle = wrapAngle(angle)
    angle = atan2(sin(angle), cos(angle));
end
