% est_traj.csv(또는 동일 형식 파일)에 저장된 pose 로그로부터
% 각 시점 사이의 상대 odometry를 계산한다.
%
% 좌표계/표기 규약:
%   p_l(k)   : k시점 위치. OpenVINS local frame에서 표현된 값
%   R_l_i(k) : k시점 IMU 자세. OpenVINS local frame에 대한 IMU/body 자세
%
% 따라서 상대 이동량은 다음과 같이 계산된다.
%   dp_local(k) = p_l(k) - p_l(k-1)
%   dp_body(k)  = R_l_i_yaw(k-1)' * dp_local(k)
%   dyaw(k)     = yaw(k) - yaw(k-1)
%
% Particle filter의 state는 mapLocal/GT yaw0 local frame에 두고,
% control input은 매 step body_{k-1} frame 기준으로 넘긴다.

clear; clc;

% 데이터 불러오기
csv_path = 'D:\Comple Urban\Urban26\est_traj.ods';
opts = detectImportOptions(csv_path, 'VariableNamingRule', 'preserve');
T = readtable(csv_path, opts);

timestamp = T{:, 1};
pos_l = T{:, 2:4};      % 지역 좌표계 위치 [tx ty tz]
quat_xyzw = T{:, 5:8};  % ROS 형식 쿼터니언 [qx qy qz qw]
quat_wxyz = [quat_xyzw(:, 4), quat_xyzw(:, 1:3)]; % 내부 계산은 [w x y z] 순서를 사용한다.
quat_wxyz = normalize_quat_rows(quat_wxyz);

% 상대 오도메트리 계산
N = size(pos_l, 1);

% Delta time 계산
rel_timestamp_prev = timestamp(1:end-1);
rel_timestamp_curr = timestamp(2:end);
dt = rel_timestamp_curr - rel_timestamp_prev;

dp_local = zeros(N-1, 3);
dp_body = zeros(N-1, 2);
dquat_wxyz = zeros(N-1, 4);
eul_zyx = zeros(N-1, 3); % [yaw pitch roll]
dyaw_se2 = zeros(N-1, 1);

for k = 2:N
    p_prev = pos_l(k-1, :).';
    p_curr = pos_l(k, :).';

    q_prev = quat_wxyz(k-1, :);
    q_curr = quat_wxyz(k, :);

    R_prev = quat_wxyz_to_rotm(q_prev);
    R_curr = quat_wxyz_to_rotm(q_curr);

    dp_l = p_curr - p_prev;
    dR = R_prev.' * R_curr;

    yaw_prev = rotm_to_yaw(R_prev);
    yaw_curr = rotm_to_yaw(R_curr);
    R_local_body_prev = [cos(yaw_prev), -sin(yaw_prev);
                         sin(yaw_prev),  cos(yaw_prev)];
    dp_b = R_local_body_prev' * dp_l(1:2);
    dyaw = atan2(sin(yaw_curr - yaw_prev), cos(yaw_curr - yaw_prev));

    dq_wxyz = rotm_to_quat_wxyz(dR);
    dq_wxyz = dq_wxyz / norm(dq_wxyz);

    dp_local(k-1, :) = dp_l';
    dp_body(k-1, :) = dp_b';
    dquat_wxyz(k-1, :) = dq_wxyz;
    eul_zyx(k-1, :) = rotm_to_eul_zyx(dR);
    dyaw_se2(k-1) = dyaw;
end

% 결과를 다시 ROS 스타일 [x y z w] 순서로 변환한다.
dquat_xyzw = [dquat_wxyz(:, 2:4), dquat_wxyz(:, 1)];

rel_odom = table;
rel_odom.idx_prev = (1:N-1).';
rel_odom.idx_curr = (2:N).';
rel_odom.t_prev = rel_timestamp_prev;
rel_odom.t_curr = rel_timestamp_curr;
rel_odom.dt = dt;

% OpenVINS local frame 기준 병진 변화량. 디버깅/비교용.
rel_odom.dtx_local = dp_local(:, 1);
rel_odom.dty_local = dp_local(:, 2);
rel_odom.dtz_local = dp_local(:, 3);

% body_{k-1} yaw frame 기준 병진 control. Particle filter propagation용.
rel_odom.dtx_body = dp_body(:, 1);
rel_odom.dty_body = dp_body(:, 2);

% 상대 회전량
rel_odom.dqx = dquat_xyzw(:, 1);
rel_odom.dqy = dquat_xyzw(:, 2);
rel_odom.dqz = dquat_xyzw(:, 3);
rel_odom.dqw = dquat_xyzw(:, 4);

% SE(2) propagation용 yaw control.
rel_odom.dyaw = dyaw_se2;

% full relative rotation(dquat)에서 추출한 디버깅용 오일러 각.
rel_odom.dyaw_full = eul_zyx(:, 1);
rel_odom.dpitch = eul_zyx(:, 2);
rel_odom.droll = eul_zyx(:, 3);

out_path = fullfile(fileparts(csv_path), 'relative_odom.csv');
writetable(rel_odom, out_path);

fprintf('Saved relative odometry to:\n%s\n\n', out_path);
fprintf('Rows in input          : %d\n', size(T, 1));
fprintf('Rows used for odometry : %d\n', N);
fprintf('Rows in output         : %d\n', height(rel_odom));

function q = normalize_quat_rows(q)
    n = sqrt(sum(q.^2, 2));
    q = q ./ n;
end

function R = quat_wxyz_to_rotm(q)
    q = q / norm(q);
    w = q(1); x = q(2); y = q(3); z = q(4);

    R = [1 - 2*(y^2 + z^2),   2*(x*y - z*w),       2*(x*z + y*w);
         2*(x*y + z*w),       1 - 2*(x^2 + z^2),   2*(y*z - x*w);
         2*(x*z - y*w),       2*(y*z + x*w),       1 - 2*(x^2 + y^2)];
end

function q = rotm_to_quat_wxyz(R)
    tr = trace(R);

    if tr > 0
        S = sqrt(tr + 1.0) * 2;
        qw = 0.25 * S;
        qx = (R(3,2) - R(2,3)) / S;
        qy = (R(1,3) - R(3,1)) / S;
        qz = (R(2,1) - R(1,2)) / S;
    elseif (R(1,1) > R(2,2)) && (R(1,1) > R(3,3))
        S = sqrt(1.0 + R(1,1) - R(2,2) - R(3,3)) * 2;
        qw = (R(3,2) - R(2,3)) / S;
        qx = 0.25 * S;
        qy = (R(1,2) + R(2,1)) / S;
        qz = (R(1,3) + R(3,1)) / S;
    elseif R(2,2) > R(3,3)
        S = sqrt(1.0 + R(2,2) - R(1,1) - R(3,3)) * 2;
        qw = (R(1,3) - R(3,1)) / S;
        qx = (R(1,2) + R(2,1)) / S;
        qy = 0.25 * S;
        qz = (R(2,3) + R(3,2)) / S;
    else
        S = sqrt(1.0 + R(3,3) - R(1,1) - R(2,2)) * 2;
        qw = (R(2,1) - R(1,2)) / S;
        qx = (R(1,3) + R(3,1)) / S;
        qy = (R(2,3) + R(3,2)) / S;
        qz = 0.25 * S;
    end

    q = [qw, qx, qy, qz];
    q = q / norm(q);

    % Keep a consistent sign convention for easier debugging/output diffs.
    if q(1) < 0
        q = -q;
    end
end

function eul = rotm_to_eul_zyx(R)
    s = -R(3,1);
    s = max(min(s, 1.0), -1.0);
    pitch = asin(s);

    if abs(s) < 1 - 1e-12
        roll = atan2(R(3,2), R(3,3));
        yaw = atan2(R(2,1), R(1,1));
    else
        % 짐벌락 근처에서는 안정적인 대체 식을 사용한다.
        roll = atan2(-R(2,3), R(2,2));
        yaw = 0.0;
    end

    eul = [yaw, pitch, roll];
end

function yaw = rotm_to_yaw(R)
    yaw = atan2(R(2,1), R(1,1));
end
