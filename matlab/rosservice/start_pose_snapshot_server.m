function server = start_pose_snapshot_server(masterURI, nodeHost)
% ROS1 MATLAB service server for OpenVINS pose snapshot requests.
% 이 함수는 OpenVINS C++ 노드가 보낸 pose snapshot 서비스 요청을
% MATLAB에서 받아서 검증한 뒤, base workspace의 pose와 meas에 누적 저장한다.
    
    %% ROS Service Configuration
    SERVICE_NAME = "/matlab/pose_snapshot";
    SERVICE_TYPE = "ov_msckf/PoseSnapshotToMatlab";

    try
        % MATLAB이 이미 ROS에 연결되어 있으면 rosinit을 다시 하지 않는다.
        rosnode list;
    catch
        % 아직 연결되어 있지 않으면 rosinit으로 ROS 노드를 시작한다.
        rosinit(masterURI, "NodeHost", nodeHost);
    end

    % Pose와 Meas는 베이스 작업공간에 명시적으로 저장한다.
    assignin("base", "pose", struct([]));
    assignin("base", "meas", struct([]));

    % ROS 서비스 서버를 생성한다. 
    server = rossvcserver(SERVICE_NAME, SERVICE_TYPE, @poseSnapshotCallback, ...
        "DataFormat", "struct");

    fprintf("[MATLAB] Service ready: %s [%s]\n", SERVICE_NAME, SERVICE_TYPE);
    fprintf("[MATLAB] Unique accepted requests will be appended to 'pose' and 'meas'.\n");

    function resp = poseSnapshotCallback(~, req, resp)
        try
            %% =========================== Request =============================
            % Request 필드 double 배열로 변환
            snapshot = struct();
            snapshot.t_cam = double(req.TimestampCam);
            snapshot.t_imu = double(req.TimestampImu);
            snapshot.position = double(req.Position(:));
            snapshot.quaternion = double(req.Quaternion(:)); % [x y z w]
            snapshot.pose_covariance_row_major = double(req.PoseCovarianceRowMajor(:));
            snapshot.pose_covariance = reshape(snapshot.pose_covariance_row_major, 6, 6).';

            % Snapshot의 필드가 비정상이면 C++ 클라이언트에 거절 응답을 보내고 종료한다.
            [ok, errMsg] = validateSnapshot(snapshot);  
            if ~ok
                resp = writeResponse(resp, computeConstraint([], false, errMsg));
                fprintf(2, "[MATLAB] Rejected snapshot: %s\n", errMsg);
                return;
            end
            
            % snapshot을 누적하여 저장
            isNew = stackSnapshot(snapshot);

            % MATLAB 콘솔에 주요 pose 정보를 출력한다.
            fprintf("  frame = %d\n", evalin("base", "numel(pose)"));
            fprintf("  t_cam = %.9f\n", snapshot.t_cam);
            % fprintf("  t_imu = %.9f\n", snapshot.t_imu);
            % fprintf("  pos   = [%.6f %.6f %.6f]\n", snapshot.position);
            % fprintf("  quat  = [%.6f %.6f %.6f %.6f]\n", snapshot.quaternion);
            % ==================================================================
            
            %% =========================== Response =============================
            accepted = true;
            constraint = computeConstraint(snapshot, accepted, "Accepted", isNew);

            % Constraint의 필드가 비정상이면 C++ 클라이언트에 거절 응답을 보내고 종료한다.
            [accepted, errMsg] = validateConstraint(constraint);
            if ~accepted
                fprintf(2, "[MATLAB] Rejected response: %s\n", errMsg);
                constraint = computeConstraint(snapshot, accepted, errMsg);
            end

            resp = writeResponse(resp, constraint);
            % ==================================================================

        catch ME
            errMsg = sprintf("MATLAB callback error: %s", ME.message);
            fprintf(2, "[MATLAB] %s\n", errMsg);
            accepted = false;
            resp = writeResponse(resp, computeConstraint([], accepted, errMsg));        
        end
    end

end

function isNew = stackSnapshot(snapshot)
    if evalin("base", "exist('pose', 'var')")
        pose = evalin("base", "pose");
    else
        pose = struct([]);
    end

    isNew = true;
    if isempty(pose)
        pose = snapshot;
    else
        if abs(pose(end).t_cam - snapshot.t_cam) <= 1e-9
            isNew = false;
            warning("ROS Server 통신 문제로 인한 중복 Request")
            return;
        end
        pose(end + 1) = snapshot;
    end

    assignin("base", "pose", pose);
end

function stackMeasurement(measurement, isNew)
    if ~isNew
        return
    end

    if evalin("base", "exist('meas', 'var')")
        meas = evalin("base", "meas");
    else
        meas = struct([]);
    end

    if isempty(meas)
        meas = measurement;
    else
        meas(end + 1) = measurement;
    end

    assignin("base", "meas", meas);
end

function constraint = computeConstraint(snapshot, accepted, statusMessage, isNew)
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
    % ========================= Validate ============================
    if nargin < 2
        accepted = false;
    end
    if nargin < 3
        if accepted
            statusMessage = "snapshot accepted";
        else
            statusMessage = "snapshot rejected";
        end
    end
    if nargin < 4
        isNew = true;
    end

    constraint = struct();
    constraint.accepted = logical(accepted);
    constraint.status_message = char(statusMessage);

    if ~constraint.accepted
        constraint.residual = 0;
        constraint.jacobian = zeros(1, 6);
        constraint.measurement_cov = 1;
        return;
    end
    % ===============================================================

    % =========================== TO-DO =============================
    
    % %%% [오프라인] Initialization
    % % 1. Configuration 불러오기
    % cfg = config();
    % 
    % % 2. 초기값(Global pose) 불러오기
    % globalPose0 = cfg.initialPose;
    % 
    % % 3. 도로중심선 지도 불러오기
    % map = loadSDMap(cfg, globalPose0);
    %
    % % 4. 초기 자세 결정
    % globalPose0 = alignment(cfg, globalPose0, map);

    %%% [온라인] Residual 및 Jacobian 계산
    % 1. t_cam에 대응하는 Bezier Inference ('timestamp'.yaml) 불러오기
    cfg.yamlDir = 'D:\Comple Urban\Urban26\stereo_left_bezier_gt';
    cfg.bevW = 120;
    cfg.bevH = 120; 
    cfg.resolution = 0.5;
    cfg.u0 = 60;
    cfg.v0 = 120; 
    cfg.curveType = 'bezier';

    yamlName = searchYaml(cfg, snapshot.t_cam);
    meas = loadYaml(cfg, yamlName);
    stackMeasurement(meas, isNew);

    % % 2. Map matching
    % [r, H, R] = mapMatching(cfg, map, meas);
    
    % 3. Output
    constraint.residual = 0;
    constraint.jacobian = zeros(1, 6);
    constraint.measurement_cov = 1;
    % ===============================================================

end

function resp = writeResponse(resp, constraint)
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

function [ok, errMsg] = validateSnapshot(snapshot)
    % ROS request에서 받은 값들이 MATLAB 분석에 사용할 수 있는 형태인지 확인한다.
    % 실패 시 ok=false와 사람이 읽을 수 있는 오류 메시지를 반환한다.
    ok = false;
    errMsg = '';

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

    ok = true;
end

function [ok, errMsg] = validateConstraint(constraint)
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

function values = rowMajorVector(matrix)
    % Convert an m x n MATLAB matrix to the row-major vector consumed by C++.
    values = reshape(double(matrix).', [], 1);
end
