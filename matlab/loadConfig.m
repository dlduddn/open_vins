function config = loadConfig()
    matlabDir = fileparts(mfilename('fullpath'));
    repoDir = fileparts(matlabDir);

    % Raw data
    config.imageDir  = 'D:\Comple Urban\Urban26\image\stereo_left';
    config.gtDir     = 'D:\Comple Urban\Urban26\global_pose.csv';
    config.estDir    = 'D:\Comple Urban\Urban26\est_odom.txt';
    config.estRelDir = 'D:\Comple Urban\Urban26\relative_odom.csv';
    config.yamlLogProgress = false;

    % SD-Map
    config.shpPath       = 'D:\Comple Urban\Urban26\GIS\Urban26.shp';
    config.shpCRS        = 'epsg5179';
    config.utmZone       = 52;
    config.ds            = 0.5; % 리샘플링 간격 (m)
    config.showRoadTypes = {'RDD000', 'RDD001', 'RDD002', 'RDD003', 'RDD008', 'RDD009'};
    config.mapInitErrorSeed = 1; % set [] to sample a new initial map error every run
    mapInitRand = mapInitErrorRandomStream(config.mapInitErrorSeed);
    config.mapInitXError   = 0 * randn(mapInitRand, 1);          % SD Map 기준 pose x 오차, global/UTM east [m]
    config.mapInitYError   = 0 * randn(mapInitRand, 1);          % SD Map 기준 pose y 오차, global/UTM north [m]
    config.mapInitYawError = deg2rad(0) * randn(mapInitRand, 1); % SD Map 기준 pose yaw 오차 [rad]

    % IPM
    config.yamlDir        = 'D:\Comple Urban\Urban26\stereo_left_bezier_inference';
    config.curveType      = 'bezier';
    config.bevW           = 120;
    config.bevH           = 120;
    config.resolution     = 0.5;        % m/pixel
    config.u0             = 60;
    config.v0             = 120;
    config.sampleSpacing  = 0.5;        % 샘플 간 거리 (m), 선 길이에 따라 샘플 수 자동 결정

    % Camera FOV gate for curve-map association
    config.cameraKalibrYaml = fullfile(repoDir, 'config', 'kaist', 'kalibr_imucam_chain.yaml');
    config.cameraKalibrCam  = 'cam0'; % stereo_left. Use 'cam1' for stereo_right.
    [config.cameraIntrinsics, config.cameraResolution] = ...
        loadKalibrCameraConfig(config.cameraKalibrYaml, config.cameraKalibrCam);
    config.cameraYawOffsetDeg = 0.0;      % positive: camera optical axis points left of body x
    config.cameraFovMarginDeg = 2.0;      % angular slack around intrinsic FOV
    config.cameraFovMinRange  = 0.5;      % minimum body-plane range [m]
    config.cameraFovMaxRange  = config.v0 * config.resolution; % metric range is not defined by K alone
    config.cameraFovGateMode  = 'nearest'; % 'nearest' fast, 'candidate' exact FOV-limited map search
    
    % Likelihood-based initialization
    config.initUse = false;                        % false: x0 = [0; 0; 0]
    config.initUseLongitudinalCorrection = false; % false: x0 = [0; estY; estYaw]
    config.initUseTrajectoryAlignment = false;    % false: use a single frame
    config.initTrajectoryWindowFrames = 450;
    config.initTrajectoryMaxFrames = 12;

    config.initMapCropMargin    = 20.0; % local map crop margin [m]

    config.initSearchX         = 10.0;  % x search half-width [m]
    config.initSearchY         = 10.0;  % y search half-width [m]
    config.initSearchYawDeg    = 30.0;  % yaw search half-width [deg]
    config.initCoarseStepXY     = 1.0;  % coarse grid xy step [m]
    config.initCoarseStepYawDeg = 2.0;  % coarse grid yaw step [deg]

    config.initUseFineSearch    = true;
    config.initFineSerachXY     = 1.0;  % fine search half-width around coarse best [m]
    config.initFineSearchYawDeg = 2.0;  % fine search half-width around coarse best [deg]
    config.initFineStepXY       = 0.25; % fine grid xy step [m]
    config.initFineStepYawDeg   = 0.5;  % fine grid yaw step [deg]

    config.initUsePriorScore = false;
    config.initPriorStdX = 3.0;       % weak prior for poorly observable forward shift [m]
    config.initPriorStdY = 6.0;       % lateral shift prior [m]
    config.initPriorStdYawDeg = 15.0; % yaw prior [deg]
    config.initBatchSize = 512;

    % Particle filter
    config.pfRandomSeed        = 1;
    config.pfNumParticles      = 500;
    
    % P0
    config.pfUseMeasurementUpdate = true; % false: propagation only

    config.pfInitXStd          = 3;   % [m] 3
    config.pfInitYStd          = 1.5; % [m] 1.5
    config.pfInitYawStdDeg     = 2;   % [deg] 2
    
    config.pfResampleRatio     = 0.5;
    config.pfLikelihoodTemperature = 6.0; % Tempering & Clipping
    config.pfLikelihoodMaxLogSpan  = 22.0;
    config.pfWeightUniformMix      = 0.02; 
    config.pfMinAssociatedParticleRatio = 0.02;

    % Ablation-selected combo: map-frame process noise + no query-only filters + no body-road gate.
    config.pfProcessNoiseFrame     = 'body';
    config.pfProcessForwardStd     = 0.04;  % per-step body x noise floor [m]
    config.pfProcessLateralStd     = 0.04;  % per-step body y noise floor [m]
    config.pfProcessYawStdDeg      = 0.10;  % per-step yaw noise floor [deg]
    config.pfProcessTransScale     = 0.015; % extra xy noise per meter traveled
    config.pfProcessYawScale       = 0.03;  % extra yaw noise per abs yaw input

    config.pfRoughenAfterResample  = false;
    config.pfRoughenXStd           = 0.01; % [m]
    config.pfRoughenYStd           = 0.02; % [m]
    config.pfRoughenYawStdDeg      = 0.03; % [deg]
    
    config.pfLogProgress           = true;
    config.pfLogInterval           = 50;

    config.pfVisualizeMapMatching = false; % show real-time PF curve-map matching
    config.pfVisualizeInterval = 10;       % draw every N frames
    config.pfVisualizeMapWindow = 45.0;   % half-width around estimate [m]
    config.pfVisualizeMaxParticles = 300; % cap plotted particles for speed

    % Bezier curve-map likelihood / data association
    config.assocCurveSampleSpacing = 0.5;  % Query Bezier curve sampling interval [m]
    config.assocMinCurveSamples    = 0;    % minimum samples per cubic curve
    config.assocMaxCurveSamples    = inf;  % maximum samples per cubic curve

    config.assocUseYawingRejection    = false;
    config.assocUsePitchingRejection  = false;
    config.assocMaxYawRateDegPerSec   = 25.0; % yawing frame rejection threshold [deg/s]
    config.assocMaxPitchRateDegPerSec = 10.0; % pitching frame rejection threshold [deg/s]

    config.assocUseMinCurveChordRule   = true;
    config.assocMinCurveChord          = 4.0;  % reject if curve start/end are too close [m]

    config.assocUseMinBodyDistanceRule = false;
    config.assocMinBodyDistance        = 4.0;  % reject if curve end point gets too close to body origin [m]

    config.assocUseStartConnectionRule = false;
    config.assocRequireStartConnection = true; % single-curve frames are valid lane observations
    config.assocAllowSingleCurveWithoutConnection = true;
    config.assocConnectionTol      = 3.0;  % curve endpoint connection tolerance [m]

    config.assocUseBodyRoadGate    = false; % reject particles outside road-width gate
    config.assocBodyGateBack       = 2.0;  % map points behind body allowed for road gate [m]
    config.assocBodyGateLookahead  = 30.0; % map points ahead used for road gate [m]

    config.assocUseCurveMapCorrespondenceRule = true;
    config.assocDefaultRoadWidth   = 6.0;  % fallback road width if RVWD is missing [m]
    config.assocRoadWidthScale     = 1.5;  % body-map lateral gate scale
    config.assocMinCurveMatchFraction = 0.45; % minimum matched samples for one curve
    
    config.assocMaxCurveMapDist    = 20.0;  % max distance for curve-map correspondence [m]
    config.assocMinMatchedCurves   = 1;    % minimum map-associated query curves

    config.assocUseMapKdTree       = true; % use KD-tree nearest-neighbor acceleration
    config.assocUseCameraFovGate   = true; % match only map centerlines inside camera FOV

    config.likelihoodUseRobust     = false;
    config.likelihoodSigma         = 1.75; % centerline residual std [m]
    config.likelihoodRobustScale   = 1.0;  % Cauchy-style residual scale multiplier
    config.likelihoodMaxEffectiveSamples = inf; % cap per-curve information
    config.likelihoodCurveMissLogPenalty = -24.0; % per usable curve when no association exists
    config.likelihoodMaxCurvePenalty = 20.0; % accepted matches are never worse than this
    config.likelihoodInlierReward  = 0.75; % small reward for high inlier fraction
    config.likelihoodMissLogPenalty = -90.0; % candidate penalty when no association exists

    % Ablation experiment bookkeeping
    config.ablationTag = "";
    config.ablationNotes = "";
end

function stream = mapInitErrorRandomStream(seed)
    if isempty(seed)
        stream = RandStream.getGlobalStream();
    else
        stream = RandStream('mt19937ar', 'Seed', seed);
    end
end
