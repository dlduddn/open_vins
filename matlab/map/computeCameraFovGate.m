function gate = computeCameraFovGate(cfg)
%COMPUTECAMERAFOVGATE Build a body-plane FOV gate from camera intrinsics.
%
% Camera intrinsics define angular FOV. The metric near/far range must come
% from BEV limits or explicit configuration because K alone has no scale.

    gate = struct( ...
        'enabled', false, ...
        'leftAngle', NaN, ...
        'rightAngle', NaN, ...
        'yawOffset', 0.0, ...
        'margin', 0.0, ...
        'minRange', 0.0, ...
        'maxRange', inf, ...
        'mode', 'nearest');

    if ~getLogical(cfg, 'assocUseCameraFovGate', false)
        return;
    end

    if ~isfield(cfg, 'cameraIntrinsics') || numel(cfg.cameraIntrinsics) < 4 || ...
            ~isfield(cfg, 'cameraResolution') || numel(cfg.cameraResolution) < 2
        return;
    end

    intr = double(cfg.cameraIntrinsics(:).');
    resolution = double(cfg.cameraResolution(:).');
    fx = intr(1);
    cx = intr(3);
    imageW = resolution(1);

    if ~isfinite(fx) || fx <= 0 || ~isfinite(cx) || ~isfinite(imageW) || imageW <= 0
        return;
    end

    gate.leftAngle = atan(max(cx, 0) / fx);
    gate.rightAngle = atan(max(imageW - cx, 0) / fx);
    gate.yawOffset = deg2rad(getScalar(cfg, 'cameraYawOffsetDeg', 0.0));
    gate.margin = deg2rad(getScalar(cfg, 'cameraFovMarginDeg', 0.0));
    gate.minRange = getScalar(cfg, 'cameraFovMinRange', 0.0);

    defaultMaxRange = inf;
    if isfield(cfg, 'v0') && isfield(cfg, 'resolution')
        defaultMaxRange = cfg.v0 * cfg.resolution;
    end
    gate.maxRange = getScalar(cfg, 'cameraFovMaxRange', defaultMaxRange);
    gate.mode = lower(string(getField(cfg, 'cameraFovGateMode', 'nearest')));
    gate.enabled = isfinite(gate.leftAngle) && isfinite(gate.rightAngle) && ...
        gate.leftAngle > 0 && gate.rightAngle > 0 && gate.maxRange > gate.minRange;
end
