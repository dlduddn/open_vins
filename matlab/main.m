%% Init
clear all
close all
clc

addpath estimator\ map\ math\ rosservice\

%% Configuration
config = loadConfig();

%% Preprocessing_1: VIO Estimate & Ground Truth
% Estimated SE(2) in 'OpenVINS local' frame
estRawData = readtable(config.estDir);
epoch   = estRawData{:,1}';

est.Pos  = [estRawData{:,2}, estRawData{:,3}]';
est.Pos0  = [estRawData{1,2}, estRawData{1,3}]';
est.Quat0 = [estRawData{1,8}, estRawData{1,5}, estRawData{1,6}, estRawData{1,7}]'; % [qw qx qy qz]
est.RPY0  = Quat2RPY(est.Quat0);
est.Yaw0  = est.RPY0(3);

est.R_ovlocal_body0 = [cos(est.Yaw0), -sin(est.Yaw0);
                   sin(est.Yaw0),  cos(est.Yaw0)];
est.PoseBody0 = est.R_ovlocal_body0' * (est.Pos - est.Pos0); % (Note: OpenVINS의 local은 Body0가 아님)

% Relative SE(2) odometry in body_{k-1} frame, to be used for propagation in particle filter
estRelRawData = readtable(config.estRelDir);
len = height(estRelRawData) + 1;
estRel.dYaw  = atan2(sin(estRelRawData.dyaw), cos(estRelRawData.dyaw))';
estRel.dSE2 = [estRelRawData.dtx_body'; estRelRawData.dty_body'; estRel.dYaw]; % 3 x (len-1)

% True SE(3) in 'global/UTM' frame
gt.rawData = readmatrix(config.gtDir); % [r00 r01 r02 tx r10 r11 r12 ty r20 r21 r22 tz]
gt.Timestamps  = gt.rawData(:, 1);
gt.PoseMapGlobal = pose12ToSE3(gt.rawData(:, 2:13));

gt.PoseMapLocal = transformToInitLocal(gt.PoseMapGlobal, "Yaw");
gt.PosMapLocal = reshape(gt.PoseMapLocal(1:2, 4, :), 2, []);

%% Preprocessing_2: Bezier Curve
% epoch(k)에 대응하는 VIO image Bezier inference를 불러와 measurement로 저장한다.
bezierQuery = loadBezier(config, epoch); % meas.pointBody = [xCam yCam width id]

%% Preprocessing_3: Road Centerline Map in 'map local' frame
mapDB = buildMap(config, gt.PoseMapGlobal(:, :, 1)); % [xMap_Local yMap_Local 차선수 도로폭]

%% Initialization_1: Map alignment
% buildMap을 좌표계 변환하는 gt.PoseMapGlobal(:, :, 1)에 오차가 첨가 되지 않으면 x0 = [0 0 0] 이다.
% 하지만, 실제상황에서는 절대 참 위치는 모르고 대략적인 초기값이 주어진다. 이 상황에서 meas와 map 간의 정렬을 통해 
% x0을 찾는다

%% Initialization_2: State
% Initial VIO pose error in mapLocal frame: [x; y; yaw].
% This can be replaced by the result of Initialization_1 map alignment.
trueInitStateMapLocal = [ ...
    gt.PosMapLocal(:,1); ...
    atan2(gt.PoseMapLocal(2,1,1), gt.PoseMapLocal(1,1,1)) ...
]; % evaluation/debug only

initXStd = 0.001;              % [m]
initYStd = 0.001;              % [m]
initYawStd = deg2rad(0.0001);  % [rad]
P0 = diag([initXStd^2, initYStd^2, initYawStd^2]);
x0 = [0; 0; 0];                % prior mean; replace with map-alignment result

%% Particle filtering
N = 1000; % particle
sqrtQ = []; % sqrtQ = chol(Q, 'lower');
R = [];

[xhat, N_eff] = sir(x0, P0, estRel.dSE2, bezierQuery, sqrtQ, R, len, N);

%%
figure;
scatter(mapDB(:,1), mapDB(:,2), 1);
hold on;

plot(gt.PosMapLocal(1,:), gt.PosMapLocal(2,:), 'k-', 'LineWidth', 1.0);
plot(est.PoseBody0(1,:), est.PoseBody0(2,:), 'b-', 'LineWidth', 1.0);
plot(xhat(1,:), xhat(2,:), 'r-', 'LineWidth', 1.0 )

scatter(gt.PosMapLocal(1,1),   gt.PosMapLocal(2,1),   30, 'g', 'filled');
scatter(gt.PosMapLocal(1,end), gt.PosMapLocal(2,end), 30, 'y', 'filled');

legend("Map", "GT", "Ref (VIO)", "Proposed", "GT Start", "GT End");

axis equal;
grid on;
xlabel("x mapLocal [m]");
ylabel("y mapLocal [m]");
title("SD Map and Trajectory in mapLocal Frame");

%% Evalusation
% % NaN값 처리
% squared_errors_SIS = (xhat_SIS - trueX).^2;
% sum_valid_squared_errors_SIS = sum(squared_errors_SIS, 3, 'omitnan');
% valid_counts_SIS = sum(~isnan(squared_errors_SIS), 3);
% xrmseSIS = sum_valid_squared_errors_SIS ./ valid_counts_SIS;
%
% squared_errors_SIR = (xhat - trueX).^2;
% sum_valid_squared_errors_SIR = sum(squared_errors_SIR, 3, 'omitnan');
% valid_counts_SIR = sum(~isnan(squared_errors_SIR), 3);
% xrmseSIR = sum_valid_squared_errors_SIR ./ valid_counts_SIR;
%
% % % NaN값 미처리
% % xrmseSIS = sum((xhat_SIS - trueX).^2, 3)/MC;
% % xrmseSIR = sum((xhat_SIR - trueX).^2, 3)/MC;
%
% figure
% subplot(2,2,1)
% plot(xrmseSIS(1,:));
% hold on;
% plot(xrmseSIR(1,:));
% title("Rmse : Position 'x' ", 'FontSize', 14)
% legend("SIS","SIR")
%
% subplot(2,2,2)
%
% plot(xrmseSIS(2,:));
% hold on;
% plot(xrmseSIR(2,:));
% title("Rmse : Position 'y'", 'FontSize', 14)
% legend("SIS","SIR")
%
% subplot(2,2,3)
%
% plot(xrmseSIS(3,:));
% hold on;
% plot(xrmseSIR(3,:));
% title("Rmse : Velocity of 'x'", 'FontSize', 14)
% legend("SIS","SIR")
%
% subplot(2,2,4)
%
% plot(xrmseSIS(4,:));
% hold on;
% plot(xrmseSIR(4,:));
% title("Rmse : Velocity of 'y'", 'FontSize', 14)
% legend("SIS","SIR")
%
% sgtitle('RMSE')
%
% %% Trajectory
% k = 1;
% figure
% plot(trueX(1,:),trueX(2,:));
% hold on
% plot(xhat_SIS(1,:,k),xhat_SIS(2,:,k), color='g');
% plot(xhat(1,:,k),xhat(2,:,k), color='r');
% title("Trajectory : True vs SIS vs SIR,     " + k+"th MC" , 'FontSize', 14)
% legend("True","SIS","SIR")
%
% %% N_eff
% figure
% plot(N_eff)
% title("Effective Sample Size", 'FontSize', 14)
