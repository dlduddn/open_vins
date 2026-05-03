function [xhat, N_eff, pfDiag] = sirLikelihood(cfg, mapDB, x0, P0, dSE2, meas, sqrtQ, len, N)
%SIRLIKELIHOOD Toy SIR particle filter using ordinary likelihood.
%
% This file intentionally does not call curveMapLikelihood().
% Measurement weights are updated directly as:
%
%     w_new = w_old .* likelihood
%
% where likelihood is computed by curveMapLikelihoodBasic().

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

    [likelihood, info] = curveMapLikelihoodBasic(cfg, mapDB, frame, xi);
    likelihood(~isfinite(likelihood) | likelihood < 0) = 0;

    matchViz.likelihood = likelihood;
    matchViz.safeLikelihood = likelihood;
    matchViz.likelihoodInfo = info;

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
