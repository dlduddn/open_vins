function config = loadConfig()
    % Raw data
    config.imageDir  = 'D:\Comple Urban\Urban26\image\stereo_left';
    config.gtDir     = 'D:\Comple Urban\Urban26\global_pose.csv';
    config.estDir    = 'D:\Comple Urban\Urban26\est_odom.txt';
    config.estRelDir = 'D:\Comple Urban\Urban26\relative_odom.csv';
    % VIO drift
    config.seed       = 1;                    % 재현성
    config.transSigma = 0;                    % 각 프레임 translation noise std [m]
    config.rotSigma   = 0;                    % 각 프레임 rotation noise std [rad]
    config.useKnownInitialPose = true;        % True: 초기자세를 알고있다고 가정하여, 오차 없음
    
    % SD-Map
    config.shpPath       = 'D:\Comple Urban\Urban26\GIS\Urban26.shp';
    config.shpCRS        = 'epsg5179';
    config.utmZone       = 52;
    config.ds            = 0.5;            % 리샘플링 간격 (m)
    config.showRoadTypes = {'RDD000', 'RDD001', 'RDD002', 'RDD003', 'RDD008', 'RDD009'};
    config.mapInitXError   = 0.0;          % SD Map 기준 pose x 오차, global/UTM east [m]
    config.mapInitYError   = 0.0;          % SD Map 기준 pose y 오차, global/UTM north [m]
    config.mapInitYawError = deg2rad(0.0); % SD Map 기준 pose yaw 오차 [rad]
    
    % IPM
    config.yamlDir        = 'D:\Comple Urban\Urban26\stereo_left_bezier_gt';
    config.bevW           = 120;
    config.bevH           = 120;
    config.resolution     = 0.5;        % m/pixel
    config.u0             = 60;
    config.v0             = 120;
    config.curveType      = 'bezier';
    config.sampleSpacing  = 0.5;        % 샘플 간 거리 (m), 선 길이에 따라 샘플 수 자동 결정
    
    % ICP
    config.windowSize            = 1;      % 현재 프레임 Bezier point만 사용
    config.icpInterval           = 1;      % ICP 수행 간격 (프레임)
    config.propagateCorrection   = false;  % Toy example: ICP 보정 결과를 다음 시점에 전파하지 않음
    config.maxIter               = 1;      % G-ICP 최대 반복
    config.tolerance             = 1e-6;   % 수렴 임계값
    config.kNeighbors            = 10;     % 공분산 추정 이웃 수
    config.covEpsilon            = 0.001;  % G-ICP 공분산 최소 고유값
    config.maxCorrespondenceDist = 15;     % 대응점 최대 거리 (m)
    config.sdCropRadius          = 150;    % 현재 위치 기준 실제 SD Map 사용 반경 (m)
    config.sourceCropMargin      = 20.0;   % source bbox 주변 SD Map 허용 여유 (m)
    config.minPoints             = 50;     % ICP 수행 최소 점 수
    config.minCorrespondences    = 30;     % 유효 대응점 최소 개수
    config.maxInitialResidual    = 10.0;   % ICP 전 source-target 중앙 최근접거리 상한
    config.maxFinalResidual      = 2.5;    % ICP 후 중앙 최근접거리 상한
    config.maxTranslationCorrection = 4.0; % 프레임별 최대 허용 translation 보정량 (m)
    config.maxRotationCorrectionDeg = 12.0; % 프레임별 최대 허용 yaw 보정량 (deg)
    config.visualize             = true;   % 매 프레임 시각화 ON/OFF
    config.vizPause              = 0;      % 프레임 간 대기 시간 (sec)
end
