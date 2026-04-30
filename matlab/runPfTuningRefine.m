function results = runPfTuningRefine(nParticles)
%RUNPFTUNINGREFINE Focused sweep around the best inference-data settings.

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

    data = loadRefineData(baseCfg);
    variants = makeRefineVariants();

    nVariant = numel(variants);
    name = strings(nVariant, 1);
    vioRmse = zeros(nVariant, 1);
    pfRmse = zeros(nVariant, 1);
    improvementPct = zeros(nVariant, 1);
    meanNeffRatio = zeros(nVariant, 1);
    elapsedSec = zeros(nVariant, 1);

    for i = 1:nVariant
        cfg = applyVariant(baseCfg, variants(i));
        fprintf('\n[%d/%d] %s: %s\n', i, nVariant, variants(i).name, variants(i).note);

        tic;
        [xhat, nEff] = sir(cfg, data.mapDB, [0; 0; 0], zeros(3), ...
            data.dSE2, data.queryBezier, [], data.len, cfg.pfNumParticles);
        elapsedSec(i) = toc;

        metrics = evaluateTrajectories(data.gtXY, data.vioXY, xhat(1:2, :));

        name(i) = variants(i).name;
        vioRmse(i) = metrics.vio.rmse;
        pfRmse(i) = metrics.pf.rmse;
        improvementPct(i) = metrics.rmseImprovementPct;
        meanNeffRatio(i) = mean(nEff, 'omitnan') / cfg.pfNumParticles;

        fprintf('  PF RMSE=%.3f m, VIO RMSE=%.3f m, improvement=%.2f%%, mean N_eff=%.2f, %.1f sec\n', ...
            pfRmse(i), vioRmse(i), improvementPct(i), meanNeffRatio(i), elapsedSec(i));
    end

    results = table(name, vioRmse, pfRmse, improvementPct, meanNeffRatio, elapsedSec);
    results = sortrows(results, 'pfRmse');

    outPath = fullfile(fileparts(mfilename('fullpath')), 'pf_tuning_refine_results.csv');
    writetable(results, outPath);
    fprintf('\nSaved refine tuning table: %s\n', outPath);
    disp(results);
end

function variants = makeRefineVariants()
    variants = struct('name', {}, 'note', {}, 'params', {});

    variants(end + 1) = variant("smooth", ...
        "stage-1 best smoother PF", smoothParams());

    variants(end + 1) = variant("smooth_low_info", ...
        "smooth PF + lower measurement information", mergeStructs(smoothParams(), struct( ...
            'pfLikelihoodTemperature', 7.0, ...
            'pfLikelihoodMaxLogSpan', 18.0, ...
            'likelihoodMaxEffectiveSamples', 10, ...
            'likelihoodInlierReward', 0.5)));

    variants(end + 1) = variant("smooth_mild_info", ...
        "smooth PF + mildly lower measurement information", mergeStructs(smoothParams(), struct( ...
            'pfLikelihoodTemperature', 6.0, ...
            'pfLikelihoodMaxLogSpan', 22.0, ...
            'likelihoodMaxEffectiveSamples', 12, ...
            'likelihoodInlierReward', 0.75)));

    variants(end + 1) = variant("smooth_very_low_info", ...
        "smooth PF + very low measurement information", mergeStructs(smoothParams(), struct( ...
            'pfLikelihoodTemperature', 9.0, ...
            'pfLikelihoodMaxLogSpan', 14.0, ...
            'likelihoodMaxEffectiveSamples', 8, ...
            'likelihoodInlierReward', 0.25)));

    variants(end + 1) = variant("smooth_strict_match", ...
        "smooth PF + require higher inlier fraction", mergeStructs(smoothParams(), struct( ...
            'assocMinCurveMatchFraction', 0.60)));

    variants(end + 1) = variant("smooth_low_info_strict_match", ...
        "smooth low-info + higher inlier fraction", mergeStructs(smoothParams(), struct( ...
            'pfLikelihoodTemperature', 7.0, ...
            'pfLikelihoodMaxLogSpan', 18.0, ...
            'likelihoodMaxEffectiveSamples', 10, ...
            'likelihoodInlierReward', 0.5, ...
            'assocMinCurveMatchFraction', 0.60)));

    variants(end + 1) = variant("smooth_strict_residual", ...
        "smooth PF + tighter residual gate", mergeStructs(smoothParams(), struct( ...
            'assocMaxCurveMapDist', 4.0, ...
            'likelihoodSigma', 1.25, ...
            'likelihoodMaxCurvePenalty', 20.0, ...
            'likelihoodCurveMissLogPenalty', -24.0)));

    variants(end + 1) = variant("smooth_soft_residual", ...
        "smooth PF + softer residual model", mergeStructs(smoothParams(), struct( ...
            'assocMaxCurveMapDist', 6.0, ...
            'likelihoodSigma', 1.75, ...
            'likelihoodMaxCurvePenalty', 20.0, ...
            'likelihoodCurveMissLogPenalty', -24.0)));

    variants(end + 1) = variant("ultra_smooth_mild_info", ...
        "even smoother PF + mild measurement information", mergeStructs(struct( ...
            'pfProcessForwardStd', 0.010, ...
            'pfProcessLateralStd', 0.030, ...
            'pfProcessYawStdDeg', 0.08, ...
            'pfProcessTransScale', 0.010, ...
            'pfProcessYawScale', 0.02, ...
            'pfRoughenXStd', 0.008, ...
            'pfRoughenYStd', 0.015, ...
            'pfRoughenYawStdDeg', 0.02), struct( ...
            'pfLikelihoodTemperature', 6.0, ...
            'pfLikelihoodMaxLogSpan', 22.0, ...
            'likelihoodMaxEffectiveSamples', 12, ...
            'likelihoodInlierReward', 0.75)));
end

function params = smoothParams()
    params = struct( ...
        'pfProcessForwardStd', 0.015, ...
        'pfProcessLateralStd', 0.04, ...
        'pfProcessYawStdDeg', 0.10, ...
        'pfProcessTransScale', 0.015, ...
        'pfProcessYawScale', 0.03, ...
        'pfRoughenXStd', 0.01, ...
        'pfRoughenYStd', 0.02, ...
        'pfRoughenYawStdDeg', 0.03);
end

function v = variant(name, note, params)
    v = struct('name', string(name), 'note', string(note), 'params', params);
end

function out = mergeStructs(varargin)
    out = struct();
    for i = 1:nargin
        names = fieldnames(varargin{i});
        for j = 1:numel(names)
            out.(names{j}) = varargin{i}.(names{j});
        end
    end
end

function cfg = applyVariant(cfg, variant)
    names = fieldnames(variant.params);
    for i = 1:numel(names)
        cfg.(names{i}) = variant.params.(names{i});
    end
end

function data = loadRefineData(cfg)
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
