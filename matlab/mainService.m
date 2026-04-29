%% One-time execution
% repoRoot = 'C:\Users\SAI_001\SynologyDrive\Code\Vision\open_vins';
% rosgenmsg(repoRoot);
% addpath(fullfile(repoRoot, 'matlab_msg_gen_ros1', 'win64', 'install', 'm'));
% savepath;
% rehash toolboxcache;
% rosmsg list

%% ROS Service Server
% MATLAB   : Server
% OpenVINS : Client
addpath mapMatching\
server = start_pose_snapshot_server("http://192.168.1.56:11311", "192.168.1.254");

%% Termination
clear all
rosshutdown

%% Results
positions = reshape([pose.position], 3, []).';
quaternions = reshape([pose.quaternion], 4, []).';
plot3(positions(:,1), positions(:,2), positions(:,3))
axis equal