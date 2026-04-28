function server = start_pose_snapshot_server(masterURI, nodeHost)
% ROS1 MATLAB service server for OpenVINS pose snapshot requests.
%
% 이 함수는 OpenVINS C++ 노드가 보낸 pose snapshot 서비스 요청을
% MATLAB에서 받아서 검증한 뒤, base workspace의 last_pose와 pose_history에 저장한다.
    
    % ROS Service Configuration
    SERVICE_NAME = "/matlab/pose_snapshot";
    SERVICE_TYPE = "ov_msckf/PoseSnapshotToMatlab";

    % MATLAB이 이미 ROS에 연결되어 있으면 rosinit을 다시 하지 않는다.
    % 아직 연결되어 있지 않으면 rosinit으로 ROS 노드를 시작한다.
    try
        rosnode list;
    catch
        rosinit(masterURI, "NodeHost", nodeHost);
    end

    % ROS callback은 MATLAB/ROS 내부 컨텍스트에서 실행되므로 상태는 base
    % workspace에 명시적으로 저장한다.
    assignin("base", "pose_history", struct([]));

    % ROS 서비스 서버를 생성한다. DataFormat을 struct로 맞춰 req/resp 필드를
    % MATLAB 구조체처럼 접근한다.
    server = rossvcserver(SERVICE_NAME, SERVICE_TYPE, @poseSnapshotCallback, ...
        "DataFormat", "struct");

    fprintf("[MATLAB] Service ready: %s [%s]\n", SERVICE_NAME, SERVICE_TYPE);
    fprintf("[MATLAB] Latest request will be stored in base workspace as 'last_pose'.\n");
    fprintf("[MATLAB] All accepted requests will be appended to 'pose_history'.\n");

    function resp = poseSnapshotCallback(~, req, resp)
        try
            % ROS service request 필드를 MATLAB에서 다루기 쉬운 double 배열로 변환한다.
            % Quaternion은 OpenVINS/JPL 관례의 [x y z w] 순서를 그대로 보존한다.
            snapshot = struct();
            snapshot.t_cam = double(req.TimestampCam);
            snapshot.t_imu = double(req.TimestampImu);
            snapshot.position = double(req.Position(:));
            snapshot.quaternion = double(req.Quaternion(:));
            snapshot.pose_covariance_row_major = double(req.PoseCovarianceRowMajor(:));
            [ok, errMsg, cov6x6] = validateSnapshot(snapshot);

            % 필수 필드가 비정상이면 C++ 클라이언트에 거절 응답을 보내고 종료한다.
            if ~ok
                resp = writeConstraintResponse(resp, emptyConstraint(false, errMsg));
                fprintf(2, "[MATLAB] Rejected snapshot: %s\n", errMsg);
                return;
            end

            % 검증이 끝난 covariance는 MATLAB 행렬 형태로 함께 저장한다.
            snapshot.pose_covariance = cov6x6;
            snapshot.received_at = datetime("now", ...
                "TimeZone", "local", ...
                "Format", "yyyy-MM-dd HH:mm:ss.SSS");

            % 사용자가 바로 확인할 수 있게 base workspace에 최신 snapshot과
            % 누적 history를 기록한다.
            history_count = appendSnapshotToHistory(snapshot);

            % MATLAB 콘솔에 주요 pose 정보를 출력한다.
            fprintf("[MATLAB] Snapshot received index=%d t_cam=%.9f t_imu=%.9f\n", ...
                history_count, snapshot.t_cam, snapshot.t_imu);
            fprintf("  pos   : [%.6f %.6f %.6f]\n", snapshot.position);
            fprintf("  quat  : [%.6f %.6f %.6f %.6f]\n", snapshot.quaternion);

            % =========================== Response =============================
            constraint = computeConstraint(snapshot);

            [ok, errMsg] = validateConstraintPayload(constraint);
            if ~ok
                fprintf(2, "[MATLAB] Rejected constraint response: %s\n", errMsg);
                constraint = emptyConstraint(false, errMsg);
            end

            resp = writeConstraintResponse(resp, constraint);
            % ==================================================================

        catch ME
            errMsg = sprintf("MATLAB callback error: %s", ME.message);
            fprintf(2, "[MATLAB] %s\n", errMsg);
            resp = writeConstraintResponse(resp, emptyConstraint(false, errMsg));
        end
    end

end

function history_count = appendSnapshotToHistory(snapshot)
    if evalin("base", "exist('pose_history', 'var')")
        pose_history = evalin("base", "pose_history");
    else
        pose_history = struct([]);
    end

    if isempty(pose_history)
        pose_history = snapshot;
    else
        pose_history(end + 1) = snapshot;
    end

    assignin("base", "last_pose", snapshot);
    assignin("base", "pose_history", pose_history);
    history_count = numel(pose_history);
end

function [ok, errMsg, cov6x6] = validateSnapshot(snapshot)
    % ROS request에서 받은 값들이 MATLAB 분석에 사용할 수 있는 형태인지 확인한다.
    % 실패 시 ok=false와 사람이 읽을 수 있는 오류 메시지를 반환한다.
    ok = false;
    errMsg = '';
    cov6x6 = [];

    % 카메라/IMU timestamp는 단일 유한 실수여야 한다.
    if ~isscalar(snapshot.t_cam) || ~isfinite(snapshot.t_cam)
        errMsg = 't_cam is invalid';
        return;
    end

    if ~isscalar(snapshot.t_imu) || ~isfinite(snapshot.t_imu)
        errMsg = 't_imu is invalid';
        return;
    end

    % 위치는 [x y z] 3개 값이어야 한다.
    if numel(snapshot.position) ~= 3 || any(~isfinite(snapshot.position))
        errMsg = 'pos must be 3 finite values';
        return;
    end

    % JPL quaternion은 [x y z w] 4개 값이어야 한다.
    if numel(snapshot.quaternion) ~= 4 || any(~isfinite(snapshot.quaternion))
        errMsg = 'quat must be 4 finite values';
        return;
    end

    % Pose covariance는 C++에서 row-major 6x6 행렬을 펼친 36개 값으로 온다.
    if numel(snapshot.pose_covariance_row_major) ~= 36 || ...
            any(~isfinite(snapshot.pose_covariance_row_major))
        errMsg = 'pose_covariance_row_major must be 36 finite values';
        return;
    end

    % C++ row-major 36-vector -> MATLAB 6x6 matrix.
    % MATLAB reshape는 column-major 기준이므로 transpose를 적용해 원래 행 순서를 복원한다.
    cov6x6 = reshape(snapshot.pose_covariance_row_major, 6, 6).';

    ok = true;
end

function constraint = computeConstraint(snapshot)
    % TODO: Replace this template with the real MATLAB-side constraint.
    %
    % C++ expects:
    %   residual: m x 1
    %   jacobian: m x 6
    %   measurement_cov: m x m
    %
    % Jacobian columns must follow:
    %   [position_error(3), orientation_error(3)]
    %
    % Residual convention should match OpenVINS EKFUpdate(), which applies:
    %   dx = K * residual
    %
    % For now this is a no-op debug constraint. It exercises the MATLAB->C++
    % response path without changing state or covariance.
    constraint = struct();
    constraint.accepted = false;
    constraint.status_message = sprintf( ...
        "snapshot received; TODO no-op constraint returned at t_cam=%.9f", ...
        snapshot.t_cam);
    constraint.residual = 0;
    constraint.jacobian = zeros(1, 6);
    constraint.measurement_cov = 1;
end

function [ok, errMsg] = validateConstraintPayload(constraint)
    % Ensure the response shape is exactly what the C++ EKF bridge expects.
    ok = false;
    errMsg = '';

    if ~isfield(constraint, "accepted") || ~isscalar(constraint.accepted)
        errMsg = 'constraint.accepted must be a scalar logical';
        return;
    end

    if ~constraint.accepted
        ok = true;
        return;
    end

    residual = double(constraint.residual(:));
    H = double(constraint.jacobian);
    R = double(constraint.measurement_cov);

    m = numel(residual);
    if m == 0
        errMsg = 'accepted constraint must have at least one residual';
        return;
    end

    if size(H, 1) ~= m || size(H, 2) ~= 6
        errMsg = 'jacobian must be m x 6';
        return;
    end

    if size(R, 1) ~= m || size(R, 2) ~= m
        errMsg = 'measurement_cov must be m x m';
        return;
    end

    if any(~isfinite(residual)) || any(~isfinite(H(:))) || any(~isfinite(R(:)))
        errMsg = 'constraint payload contains NaN or Inf';
        return;
    end

    ok = true;
end

function constraint = emptyConstraint(accepted, statusMessage)
    constraint = struct();
    constraint.accepted = logical(accepted);
    constraint.status_message = char(statusMessage);
    constraint.residual = zeros(0, 1);
    constraint.jacobian = zeros(0, 6);
    constraint.measurement_cov = zeros(0, 0);
end

function resp = writeConstraintResponse(resp, constraint)
    % MATLAB stores matrices column-major, while the ROS payload is row-major.
    residual = double(constraint.residual(:));
    H = double(constraint.jacobian);
    R = double(constraint.measurement_cov);

    resp.Accepted = logical(constraint.accepted);
    resp.StatusMessage = char(constraint.status_message);
    resp.Residual = residual;
    resp.JacobianRows = uint32(size(H, 1));
    resp.JacobianCols = uint32(size(H, 2));
    resp.JacobianRowMajor = rowMajorVector(H);
    resp.MeasurementCovRows = uint32(size(R, 1));
    resp.MeasurementCovCols = uint32(size(R, 2));
    resp.MeasurementCovRowMajor = rowMajorVector(R);
end

function values = rowMajorVector(matrix)
    % Convert an m x n MATLAB matrix to the row-major vector consumed by C++.
    values = reshape(double(matrix).', [], 1);
end
