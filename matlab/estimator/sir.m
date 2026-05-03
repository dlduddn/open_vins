function [xhat, N_eff, pfDiag] = sir(cfg, mapDB, x0, P0, dSE2, meas, sqrtQ, len, N)
%SIRBEZIERLIKELIHOOD Sequential importance resampling with Bezier likelihood.
%
% The measurement model is curveMapLikelihood(), shared with initialization.

    rng(cfg.pfRandomSeed);

    % Pre-allocation
    stateDim = 3;
    xhat = zeros(stateDim, len);
    N_eff = zeros(1, len);
    pfDiag = initPfDiagnostics(len);
    mapIndex = prepareCurveMapIndex(cfg, mapDB);

    % Initialization
    xi = repmat(x0(:), 1, N) + sqrt(P0) * randn(stateDim, N);
    xi(3, :) = wrapAngle(xi(3, :));
    wi = (1 / N) * ones(1, N);

    used = false;
    usableCount = 0;
    meanLogL = NaN;
    matchViz = emptyMatchViz();
    if cfg.pfUseMeasurementUpdate
        [wi, used, usableCount, meanLogL, matchViz] = updateWeights(cfg, mapIndex, meas, 1, xi, wi);
    end

    [xhat(:, 1), N_eff(1)] = weightedEstimate(xi, wi);

    % Log
    visualizePfMapMatching(cfg, mapIndex, 1, len, xi, wi, xhat(:, 1), N_eff(1), used, false, usableCount, meanLogL, matchViz);
    pfDiag = recordPfDiagnostics(pfDiag, 1, used, false, usableCount, meanLogL, N_eff(1), N, matchViz);
    logFrameProgress(cfg.pfLogProgress, 1, len, cfg.pfLogInterval, used, false, usableCount, meanLogL, N_eff(1), N, xhat(:, 1));

    for t = 2:len
        dX = dSE2(1, t - 1);
        dY = dSE2(2, t - 1);
        dYaw = dSE2(3, t - 1);

        c = cos(xi(3, :));
        s = sin(xi(3, :));
        xi(1, :) = xi(1, :) + c .* dX - s .* dY;
        xi(2, :) = xi(2, :) + s .* dX + c .* dY;
        xi(3, :) = wrapAngle(xi(3, :) + dYaw);

        sqrtQStep = processNoiseSqrt(cfg, sqrtQ, dX, dY, dYaw, stateDim);
        xi = addBodyFrameNoise(cfg, xi, sqrtQStep, N);

        used = false;
        usableCount = 0;
        meanLogL = NaN;
        resampled = false;
        matchViz = emptyMatchViz();

        if cfg.pfUseMeasurementUpdate
            [wi, used, usableCount, meanLogL, matchViz] = updateWeights(cfg, mapIndex, meas, t, xi, wi);
            [xhat(:, t), N_eff(t)] = weightedEstimate(xi, wi);
            resampled = N_eff(t) < cfg.pfResampleRatio * N;
            visualizePfMapMatching(cfg, mapIndex, t, len, xi, wi, xhat(:, t), ...
                N_eff(t), used, resampled, usableCount, meanLogL, matchViz);

            if resampled
                idx = sysresample(wi);
                xi = xi(:, idx);
                xi = roughenParticles(cfg, xi, N);
                wi = (1 / N) * ones(1, N);
            end
        else
            [xhat(:, t), N_eff(t)] = weightedEstimate(xi, wi);
            visualizePfMapMatching(cfg, mapIndex, t, len, xi, wi, xhat(:, t), ...
                N_eff(t), used, false, usableCount, meanLogL, matchViz);
        end

        pfDiag = recordPfDiagnostics(pfDiag, t, used, resampled, usableCount, meanLogL, N_eff(t), N, matchViz);
        logFrameProgress(cfg.pfLogProgress, t, len, cfg.pfLogInterval, used, resampled, ...
            usableCount, meanLogL, N_eff(t), N, xhat(:, t));
    end
end

function pfDiag = initPfDiagnostics(len)
    pfDiag = struct( ...
        'usedUpdate', false(1, len), ...
        'resampled', false(1, len), ...
        'usableCurves', zeros(1, len), ...
        'meanLogL', NaN(1, len), ...
        'nEff', NaN(1, len), ...
        'nEffRatio', NaN(1, len), ...
        'associatedParticleRatio', NaN(1, len), ...
        'roadGatePassRatio', NaN(1, len), ...
        'matchedCurveCountMean', NaN(1, len), ...
        'meanResidual', NaN(1, len), ...
        'rejectYawing', zeros(1, len), ...
        'rejectPitching', zeros(1, len), ...
        'rejectShortChord', zeros(1, len), ...
        'rejectNearBody', zeros(1, len), ...
        'rejectIsolatedStart', zeros(1, len));
end

function pfDiag = recordPfDiagnostics(pfDiag, t, used, resampled, usableCount, meanLogL, neff, nParticle, matchViz)
    pfDiag.usedUpdate(t) = logical(used);
    pfDiag.resampled(t) = logical(resampled);
    pfDiag.usableCurves(t) = usableCount;
    pfDiag.meanLogL(t) = meanLogL;
    pfDiag.nEff(t) = neff;
    pfDiag.nEffRatio(t) = neff / max(nParticle, 1);

    if isstruct(matchViz)
        if isfield(matchViz, 'associatedRatio')
            pfDiag.associatedParticleRatio(t) = matchViz.associatedRatio;
        end
        if isfield(matchViz, 'frame') && isstruct(matchViz.frame) && isfield(matchViz.frame, 'rejectSummary')
            rejectSummary = matchViz.frame.rejectSummary;
            pfDiag.rejectYawing(t) = getRejectCount(rejectSummary, 'yawing');
            pfDiag.rejectPitching(t) = getRejectCount(rejectSummary, 'pitching');
            pfDiag.rejectShortChord(t) = getRejectCount(rejectSummary, 'shortChord');
            pfDiag.rejectNearBody(t) = getRejectCount(rejectSummary, 'nearBody');
            pfDiag.rejectIsolatedStart(t) = getRejectCount(rejectSummary, 'isolatedStart');
        end
        if isfield(matchViz, 'likelihoodInfo') && isstruct(matchViz.likelihoodInfo)
            info = matchViz.likelihoodInfo;
            if isfield(info, 'roadGateOk') && ~isempty(info.roadGateOk)
                pfDiag.roadGatePassRatio(t) = mean(logical(info.roadGateOk));
            end
            if isfield(info, 'matchedCurveCount') && ~isempty(info.matchedCurveCount)
                pfDiag.matchedCurveCountMean(t) = mean(info.matchedCurveCount, 'omitnan');
            end
            if isfield(info, 'meanResidual') && ~isempty(info.meanResidual)
                residual = info.meanResidual(isfinite(info.meanResidual));
                if ~isempty(residual)
                    pfDiag.meanResidual(t) = mean(residual, 'omitnan');
                end
            end
        end
    end
end

function [wi, used, usableCount, meanLogL, matchViz] = updateWeights(cfg, mapDB, meas, t, xi, wi)
    % 기본 출력값 초기화
    % used가 true가 되는 경우에만 이번 측정값이 weight update에 실제 사용된 것임
    used = false;
    usableCount = 0;
    meanLogL = NaN;
    matchViz = emptyMatchViz();

    % 현재 시간 index t에 해당하는 측정값이 없으면 갱신하지 않음
    curTime = meas.Timestamp(t);
    measTime = meas.Bezier(t).t_cam;
    notSame = logical(curTime-measTime);
    if notSame
        return;
    end

    % 현재 측정값을 map matching에 사용할 Bezier curve frame으로 변환
    frame = buildQueryBezierFrame(cfg, meas.Bezier(t));
    usableCount = numel(frame.usableIdx); % 실제 likelihood 계산에 사용할 수 있는 curve 개수
    matchViz.frame = frame;
    if ~frame.hasUsableCurves % 사용할 수 있는 curve가 없으면 weight update 생략
        return; 
    end

    % 각 particle xi에 대해 현재 frame이 mapDB와 얼마나 잘 맞는지 log-likelihood 계산
    [logL, likelihoodInfo] = curveMapLikelihood(cfg, mapDB, frame, xi);
    matchViz.logL = logL;
    matchViz.likelihoodInfo = likelihoodInfo;
    if all(~isfinite(logL)) % 모든 particle의 likelihood가 NaN/Inf이면 갱신 불가
        return;
    end

    % particle들이 map curve와 얼마나 잘 association 되었는지 평가
    associatedRatio = particleAssociationRatio(cfg, likelihoodInfo);
    matchViz.associatedRatio = associatedRatio;

    % association 비율이 너무 낮으면 outlier 측정으로 보고 갱신하지 않음
    minAssociatedRatio = getScalar(cfg, 'pfMinAssociatedParticleRatio', 0.0);
    if isfinite(associatedRatio) && associatedRatio < minAssociatedRatio
        return;
    end

    % 유효한 log-likelihood만 구분
    finiteMask = isfinite(logL);
    safeLogL = logL;

    % NaN/Inf likelihood는 계산 안정성을 위해 매우 낮은 finite 값으로 대체
    % safeLogL(~finiteMask) = min(logL(finiteMask)) - 100.0;
    safeLogL(~finiteMask) = [];

    % likelihood가 너무 날카롭거나 극단적이지 않도록 tempering/clipping 적용
    % safeLogL = temperAndClipLogLikelihood(cfg, safeLogL);
    matchViz.safeLogL = safeLogL;

    % 현재 frame의 평균 log-likelihood 저장
    meanLogL = mean(safeLogL);

    % Bayesian weight update:
    % new weight ∝ old weight * likelihood
    % log domain에서는 logW = log(old weight) + log-likelihood
    logW = log(max(wi, realmin)) + safeLogL;

    % exp 계산 시 overflow를 막기 위한 수치 안정화
    logW = logW - max(logW);

    % log weight를 일반 weight로 변환
    wi = exp(logW);
    wsum = sum(wi);

    % weight 합이 비정상이면 uniform weight로 복구하고 갱신 실패로 반환
    if wsum <= 0 || ~isfinite(wsum)
        wi = (1 / numel(wi)) * ones(size(wi));
        return;
    end

    % weight 합이 1이 되도록 정규화
    wi = wi / wsum;

    % particle weight가 특정 particle에 과도하게 몰리는 것을 완화하기 위해
    % uniform distribution을 일부 섞을 수 있음
    % uniformMix = min(max(getScalar(cfg, 'pfWeightUniformMix', 0.0), 0.0), 1.0);
    % if uniformMix > 0
    %     wi = (1 - uniformMix) * wi + uniformMix / numel(wi);
    %     wi = wi / sum(wi);
    % end

    % 여기까지 도달하면 이번 측정값으로 weight update 성공
    used = true;
end

function matchViz = emptyMatchViz()
    matchViz = struct( ...
        'frame', [], ...
        'logL', [], ...
        'safeLogL', [], ...
        'likelihoodInfo', [], ...
        'associatedRatio', NaN);
end

function sqrtQStep = processNoiseSqrt(cfg, sqrtQ, dX, dY, dYaw, stateDim)
    if ~isempty(sqrtQ)
        if ~isequal(size(sqrtQ), [stateDim, stateDim])
            error('sqrtQ must be a %d x %d matrix.', stateDim, stateDim);
        end
        sqrtQStep = sqrtQ;
        return;
    end

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

function logL = temperAndClipLogLikelihood(cfg, logL)
    temperature = max(getScalar(cfg, 'pfLikelihoodTemperature', 1.0), eps);
    logL = logL / temperature;

    maxSpan = getScalar(cfg, 'pfLikelihoodMaxLogSpan', inf);
    if isfinite(maxSpan) && maxSpan > 0 && any(isfinite(logL))
        bestLogL = max(logL(isfinite(logL)));
        logL = max(logL, bestLogL - maxSpan);
    end
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

function [x, neff] = weightedEstimate(xi, wi)
    x = zeros(3, 1);
    x(1:2) = sum(xi(1:2, :) .* wi, 2);
    x(3) = atan2(sum(sin(xi(3, :)) .* wi), sum(cos(xi(3, :)) .* wi));
    neff = 1 / sum(wi .^ 2);
end

function logFrameProgress(logProgress, t, len, logInterval, used, resampled, ...
        usableCount, meanLogL, neff, nParticle, x)
    if ~logProgress || ~shouldLogFrame(t, len, logInterval)
        return;
    end

    if used
        zStatus = 'used';
    else
        zStatus = 'skip';
    end

    if resampled
        resampleStatus = 'yes';
    else
        resampleStatus = 'no';
    end

    progressPct = 100.0 * t / max(len, 1);
    fprintf(['PF Bezier [%5d/%5d %6.2f%%] z=%s, curves=%d, ' ...
             'N_eff=%.1f/%d, meanLogL=%.2f, x=[%.3f %.3f %.2fdeg], resampled=%s\n'], ...
        t, len, progressPct, zStatus, usableCount, neff, nParticle, ...
        meanLogL, x(1), x(2), rad2deg(x(3)), resampleStatus);
end

function doLog = shouldLogFrame(t, len, logInterval)
    doLog = t == 1 || t == len || mod(t, logInterval) == 0;
end
