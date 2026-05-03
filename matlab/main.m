%% Init
clear all
close all
clc

addpath estimator\ map\ math\ rosservice\ utils\

%% Configuration
config = loadConfig();

%% Preprocessing_1: VIO Estimate & Ground Truth
% Estimated SE(2) in 'OpenVINS local' frame
estRawData = readtable(config.estDir, 'VariableNamingRule', 'preserve');
query.Timestamp   = estRawData{:,1}';

est.Pos  = [estRawData{:,2}, estRawData{:,3}]';
est.Pos0  = [estRawData{1,2}, estRawData{1,3}]';
est.Quat = [estRawData{:,8}, estRawData{:,5}, estRawData{:,6}, estRawData{:,7}]'; % [qw qx qy qz]
est.Quat0 = [estRawData{1,8}, estRawData{1,5}, estRawData{1,6}, estRawData{1,7}]'; % [qw qx qy qz]
est.RPY0  = Quat2RPY(est.Quat0);
est.Yaw0  = est.RPY0(3);
est.RPY = zeros(3, size(est.Quat, 2));
for k = 1:size(est.Quat, 2)
    est.RPY(:, k) = Quat2RPY(est.Quat(:, k));
end

est.R_ovlocal_body0 = [cos(est.Yaw0), -sin(est.Yaw0);
                   sin(est.Yaw0),  cos(est.Yaw0)];
est.PoseBody0 = est.R_ovlocal_body0' * (est.Pos - est.Pos0); % (Note: OpenVINS의 local은 Body0가 아님)

% Relative SE(2) odometry in body_{k-1} frame, to be used for propagation in particle filter
estRelRawData = readtable(config.estRelDir);
len = height(estRelRawData) + 1;
estRel.dYaw  = atan2(sin(estRelRawData.dyaw), cos(estRelRawData.dyaw))';
estRel.dSE2 = [estRelRawData.dtx_body'; estRelRawData.dty_body'; estRel.dYaw]; % 3 x (len-1)

% True SE(3) in 'global/UTM' frame
[gt.Timestamps, gt.Poses] = genRef(config, query.Timestamp);
gt.PoseMapGlobal = pose12ToSE3(gt.Poses);
gt.PoseMapLocal = transformToInitLocal(gt.PoseMapGlobal, "Yaw");
gt.PosMapLocal = reshape(gt.PoseMapLocal(1:2, 4, :), 2, []);

%% Preprocessing_2: Bezier Curve
% queryTimestamp에 대응하는 Bezier inference를 불러온다.
query.Bezier = loadBezier(config, query.Timestamp); % pointsBody = [x_forward y_left width id]
query.Bezier = annotateQueryMotion(config, query.Bezier, est.RPY, query.Timestamp);

%% Preprocessing_3: Road Centerline Map in initial-pose local frame
% Initial global pose의 XY/yaw를 reference로 사용해 UTM map을 local frame으로 변환한다.
% buildMap 내부에서 cfg.mapInit*Error를 reference pose에 더해 초기 pose 오차를 모사한다.
[mapDB, mapDBGtViz] = buildMap(config, gt.PoseMapGlobal(:, :, 1)); % [xMap_Local yMap_Local 차선수 도로폭]

%% Initialization
% Query 중심선은 cubic Bezier curve로 변환한 뒤 곡선 단위 data association을 수행한다.
% 초기 pose는 Initialization 전용 임시 cost가 아니라 PF와 동일한 curve-map likelihood로 찾는다.
if getLogical(config, 'initUseTrajectoryAlignment', true)
    [x0Align, alignInfo] = estimateInitialTrajectoryAlignment(config, mapDB, query.Bezier, estRel.dSE2, len);
    x0TrueAlign = initialAlignmentTruth(config, gt.PoseMapGlobal(:, :, 1), gt.PoseMapGlobal(:, :, 1));
else
    [x0Align, alignInfo] = estimateInitialAlignment(config, mapDB, query.Bezier);
    x0TrueAlign = initialAlignmentTruth(config, gt.PoseMapGlobal(:, :, 1), gt.PoseMapGlobal(:, :, alignInfo.frameIdx));
end
visualizeInitialAlignment(mapDB, alignInfo, x0Align, x0TrueAlign);

%% Particle filtering
% Particle
a = tic;
N = config.pfNumParticles;

% Initial validation: do not inject initial pose error.
P0 = diag([config.pfInitXStd^2, config.pfInitYStd^2, deg2rad(config.pfInitYawStdDeg)^2]); % [m m rad]
x0 = selectPfInitialState(config, x0Align);

% Process noise is configured adaptively in sir() from cfg.pfProcess*.
sqrtQ = [];

[xhat, N_eff] = sirLikelihood(config, mapDB, x0, P0, estRel.dSE2, query, sqrtQ, len, N);
% [xhat, N_eff] = sir(config, mapDB, x0, P0, estRel.dSE2, query, sqrtQ, len, N);

% xhat 시각화 변환용: map frame (Noisy)에서 map frame (True)으로 가는 변환
xTrue0InMap = initialAlignmentTruth(config, gt.PoseMapGlobal(:, :, 1), gt.PoseMapGlobal(:, :, 1));
xhatGtViz = recenterStatesByPose(xhat, xTrue0InMap);
toc(a)

%% Visualization
figure;
scatter(mapDBGtViz(:,1), mapDBGtViz(:,2), 1, 'k', 'filled' );
hold on;

plot(gt.PosMapLocal(1,:), gt.PosMapLocal(2,:), 'g-', 'LineWidth', 1.0);
plot(est.PoseBody0(1,:), est.PoseBody0(2,:), 'b-', 'LineWidth', 1.0);
plot(xhatGtViz(1,:), xhatGtViz(2,:), 'r-', 'LineWidth', 1.0 )

scatter(gt.PosMapLocal(1,1),   gt.PosMapLocal(2,1),   30, 'g', 'filled');
scatter(gt.PosMapLocal(1,end), gt.PosMapLocal(2,end), 30, 'y', 'filled');

legend("Map", "GT", "Ref (VIO)", "Proposed", "GT Start", "GT End");

axis equal;
grid on;
xlabel("x GT initial local [m]");
ylabel("y GT initial local [m]");
title("SD Map and Trajectory in GT Initial Local Frame");

%% Evaluation
metrics = evaluateTrajectories(gt.PosMapLocal, est.PoseBody0, xhatGtViz(1:2, :));
fprintf('\nTrajectory position RMSE over %d samples\n', metrics.n);
fprintf('  VIO : %.3f m (axis RMSE x=%.3f, y=%.3f)\n', ...
    metrics.vio.rmse, metrics.vio.axisRmse(1), metrics.vio.axisRmse(2));
fprintf('  PF  : %.3f m (axis RMSE x=%.3f, y=%.3f)\n', ...
    metrics.pf.rmse, metrics.pf.axisRmse(1), metrics.pf.axisRmse(2));
fprintf('  Improvement: %.3f m (%.2f%%)\n', ...
    metrics.rmseImprovement, metrics.rmseImprovementPct);
if metrics.pf.rmse >= metrics.vio.rmse
    warning('PF did not beat VIO over the full trajectory. Tune likelihood/process parameters and rerun.');
end

figure
subplot(2,1,1)
plot(metrics.vio.error, 'b-', 'LineWidth', 1.0); hold on;
plot(metrics.pf.error, 'r-', 'LineWidth', 1.0);
legend("VIO", "PF");
grid on;
title("Position error norm", 'FontSize', 14)

subplot(2,1,2)
plot(N_eff, 'k-', 'LineWidth', 1.0);
grid on;
title("PF effective sample size", 'FontSize', 14)
