%% One-time execution
% repoRoot = 'C:\Users\SAI_001\SynologyDrive\Code\Vision\open_vins';
% 
% rosgenmsg(repoRoot);
% 
% addpath(fullfile(repoRoot, 'matlab_msg_gen_ros1', 'win64', 'install', 'm'));
% savepath;
% 
% clear classes;
% rehash toolboxcache;
% 
% rosmsg list

%% ROS Service Server 생성
% MATLAB: Server
% C++: Client
rosshutdown
server = start_pose_snapshot_server("http://192.168.1.56:11311", "192.168.1.254");

%% Inspect accumulated snapshots
% n = numel(pose_history);
% positions = reshape([pose_history.position], 3, []).';
% quaternions = reshape([pose_history.quaternion], 4, []).';

%% Termination
delete(server);
rosshutdown
clear all


