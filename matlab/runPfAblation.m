function results = runPfAblation(nParticles)
%RUNPFABLATION Full-trajectory ablation for the three PF convergence fixes.

    if nargin < 1 || isempty(nParticles)
        nParticles = [];
    end

    addpath estimator\ map\ math\ rosservice\

    baseCfg = loadConfig();
    baseCfg.yamlLogProgress = false;
    baseCfg.pfLogProgress = false;
    baseCfg.pfVisualizeMapMatching = false;
    if ~isempty(nParticles)
        baseCfg.pfNumParticles = nParticles;
    end

    data = loadAblationData(baseCfg);
    variants = makeVariants();

    nVariant = numel(variants);
    name = strings(nVariant, 1);
    gateFix = false(nVariant, 1);
    robustLikelihood = false(nVariant, 1);
    pfStabilization = false(nVariant, 1);
    usableFrames = zeros(nVariant, 1);
    usableCurves = zeros(nVariant, 1);
    vioRmse = zeros(nVariant, 1);
    pfRmse = zeros(nVariant, 1);
    improvement = zeros(nVariant, 1);
    improvementPct = zeros(nVariant, 1);
    meanNeffRatio = zeros(nVariant, 1);
    elapsedSec = zeros(nVariant, 1);

    for i = 1:nVariant
        cfg = configureVariant(baseCfg, variants(i));
        sqrtQ = processNoiseForVariant(variants(i));

        fprintf('\n[%d/%d] %s: gate=%d, likelihood=%d, pf=%d\n', ...
            i, nVariant, variants(i).name, variants(i).gateFix, ...
            variants(i).robustLikelihood, variants(i).pfStabilization);

        stats = countUsableBezierFrames(cfg, data.queryBezier);
        tic;
        [xhat, nEff] = sir(cfg, data.mapDB, [0; 0; 0], zeros(3), ...
            data.dSE2, data.queryBezier, sqrtQ, data.len, cfg.pfNumParticles);
        elapsedSec(i) = toc;

        metrics = evaluateTrajectories(data.gtXY, data.vioXY, xhat(1:2, :));

        name(i) = variants(i).name;
        gateFix(i) = variants(i).gateFix;
        robustLikelihood(i) = variants(i).robustLikelihood;
        pfStabilization(i) = variants(i).pfStabilization;
        usableFrames(i) = stats.usableFrames;
        usableCurves(i) = stats.usableCurves;
        vioRmse(i) = metrics.vio.rmse;
        pfRmse(i) = metrics.pf.rmse;
        improvement(i) = metrics.rmseImprovement;
        improvementPct(i) = metrics.rmseImprovementPct;
        meanNeffRatio(i) = mean(nEff, 'omitnan') / cfg.pfNumParticles;

        fprintf('  usable=%d frames/%d curves, PF RMSE=%.3f m, VIO RMSE=%.3f m, improvement=%.2f%%, %.1f sec\n', ...
            usableFrames(i), usableCurves(i), pfRmse(i), vioRmse(i), improvementPct(i), elapsedSec(i));
    end

    results = table(name, gateFix, robustLikelihood, pfStabilization, ...
        usableFrames, usableCurves, vioRmse, pfRmse, improvement, ...
        improvementPct, meanNeffRatio, elapsedSec);
    results = sortrows(results, 'pfRmse');

    outPath = fullfile(fileparts(mfilename('fullpath')), 'ablation_results.csv');
    writetable(results, outPath);
    fprintf('\nSaved ablation table: %s\n', outPath);
    disp(results);
end

function data = loadAblationData(cfg)
    estRawData = readtable(cfg.estDir, 'VariableNamingRule', 'preserve');
    queryTimestamp = estRawData{:, 1}.';

    estPos = [estRawData{:, 2}, estRawData{:, 3}].';
    estPos0 = estPos(:, 1);
    estQuat0 = [estRawData{1, 8}, estRawData{1, 5}, estRawData{1, 6}, estRawData{1, 7}].';
    estYaw0 = Quat2RPY(estQuat0);
    estYaw0 = estYaw0(3);
    r0 = [cos(estYaw0), -sin(estYaw0); sin(estYaw0), cos(estYaw0)];
    vioXY = r0.' * (estPos - estPos0);

    estQuat = [estRawData{:, 8}, estRawData{:, 5}, estRawData{:, 6}, estRawData{:, 7}].';
    estRpy = zeros(3, size(estQuat, 2));
    for k = 1:size(estQuat, 2)
        estRpy(:, k) = Quat2RPY(estQuat(:, k));
    end

    estRelRawData = readtable(cfg.estRelDir);
    len = height(estRelRawData) + 1;
    dYaw = atan2(sin(estRelRawData.dyaw), cos(estRelRawData.dyaw)).';
    dSE2 = [estRelRawData.dtx_body.'; estRelRawData.dty_body.'; dYaw];

    [~, gtPoses] = genRef(cfg, queryTimestamp);
    gtPoseMapGlobal = pose12ToSE3(gtPoses);
    gtPoseMapLocal = transformToInitLocal(gtPoseMapGlobal, "Yaw");
    gtXY = reshape(gtPoseMapLocal(1:2, 4, :), 2, []);

    queryBezier = loadBezier(cfg, queryTimestamp);
    queryBezier = annotateQueryMotion(cfg, queryBezier, estRpy, queryTimestamp);

    [mapDB, ~] = buildMap(cfg, gtPoseMapGlobal(:, :, 1));

    data = struct();
    data.len = len;
    data.dSE2 = dSE2;
    data.queryBezier = queryBezier;
    data.mapDB = mapDB;
    data.gtXY = gtXY;
    data.vioXY = vioXY;
end

function variants = makeVariants()
    variants = struct('name', {}, 'gateFix', {}, 'robustLikelihood', {}, 'pfStabilization', {});
    bools = [false, true];
    for a = bools
        for b = bools
            for c = bools
                suffix = sprintf('A%d_B%d_C%d', a, b, c);
                variants(end + 1) = struct( ...
                    'name', string(suffix), ...
                    'gateFix', a, ...
                    'robustLikelihood', b, ...
                    'pfStabilization', c); %#ok<AGROW>
            end
        end
    end
end

function cfg = configureVariant(cfg, v)
    cfg.pfRandomSeed = 7;

    if v.gateFix
        cfg.assocRequireStartConnection = false;
        cfg.assocAllowSingleCurveWithoutConnection = true;
        cfg.assocConnectionTol = 3.0;
    else
        cfg.assocRequireStartConnection = true;
        cfg.assocAllowSingleCurveWithoutConnection = false;
        cfg.assocConnectionTol = 1.5;
    end

    if v.robustLikelihood
        cfg.assocMaxCurveMapDist = 5.0;
        cfg.likelihoodUseRobust = true;
        cfg.likelihoodSigma = 1.50;
        cfg.likelihoodRobustScale = 1.0;
        cfg.likelihoodMaxEffectiveSamples = 18;
        cfg.likelihoodCurveMissLogPenalty = -28.0;
        cfg.likelihoodMaxCurvePenalty = 24.0;
        cfg.likelihoodInlierReward = 1.5;
        cfg.likelihoodMissLogPenalty = -90.0;
    else
        cfg.assocMaxCurveMapDist = 2.5;
        cfg.likelihoodUseRobust = false;
        cfg.likelihoodSigma = 0.75;
        cfg.likelihoodMaxEffectiveSamples = 25;
        cfg.likelihoodMissLogPenalty = -60.0;
    end

    if v.pfStabilization
        cfg.pfLikelihoodTemperature = 4.0;
        cfg.pfLikelihoodMaxLogSpan = 28.0;
        cfg.pfWeightUniformMix = 0.02;
        cfg.pfMinAssociatedParticleRatio = 0.02;
        cfg.pfProcessNoiseFrame = 'body';
        cfg.pfRoughenAfterResample = true;
    else
        cfg.pfLikelihoodTemperature = 1.0;
        cfg.pfLikelihoodMaxLogSpan = inf;
        cfg.pfWeightUniformMix = 0.0;
        cfg.pfMinAssociatedParticleRatio = 0.0;
        cfg.pfProcessNoiseFrame = 'map';
        cfg.pfRoughenAfterResample = false;
    end
end

function sqrtQ = processNoiseForVariant(v)
    if v.pfStabilization
        sqrtQ = [];
    else
        q = diag([3.0^2, 3.0^2, deg2rad(10.0)^2]);
        sqrtQ = chol(q, 'lower');
    end
end

function stats = countUsableBezierFrames(cfg, meas)
    usableFrames = 0;
    usableCurves = 0;
    for k = 1:numel(meas)
        frame = buildQueryBezierFrame(cfg, meas(k));
        if frame.hasUsableCurves
            usableFrames = usableFrames + 1;
            usableCurves = usableCurves + numel(frame.usableIdx);
        end
    end

    stats = struct('usableFrames', usableFrames, 'usableCurves', usableCurves);
end
