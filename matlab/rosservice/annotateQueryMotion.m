function meas = annotateQueryMotion(cfg, meas, rpy, timestamps)
%ANNOTATEQUERYMOTION Add yawing/pitching flags to Bezier measurements.
%
% Query-only data association can reject frames captured during fast yaw or
% pitch motion. This helper keeps the motion decision next to each Bezier
% frame so the likelihood model can stay measurement-centric.

    if isempty(meas)
        return;
    end

    nFrame = numel(meas);
    yawRate = zeros(1, nFrame);
    pitchRate = zeros(1, nFrame);

    if nargin >= 3 && ~isempty(rpy) && size(rpy, 1) >= 3
        yaw = rpy(3, :);
        pitch = rpy(2, :);
        if nargin < 4 || isempty(timestamps)
            timestamps = 1:nFrame;
        end
        timestamps = timestamps(:).';
        n = min([nFrame, numel(yaw), numel(pitch), numel(timestamps)]);
        timestamps = timestamps(1:n);
        dt = diff(timestamps);
        if ~isempty(dt) && median(abs(dt(isfinite(dt) & dt ~= 0))) > 1e4
            dt = dt * 1e-9;
        end
        dt = max(abs(dt), eps);

        dyaw = abs(wrapAngle(diff(yaw(1:n)))) ./ dt;
        dpitch = abs(wrapAngle(diff(pitch(1:n)))) ./ dt;

        prevYawRate = [0, dyaw];
        nextYawRate = [dyaw, 0];
        prevPitchRate = [0, dpitch];
        nextPitchRate = [dpitch, 0];

        yawRate(1:n) = max(prevYawRate, nextYawRate);
        pitchRate(1:n) = max(prevPitchRate, nextPitchRate);
    end

    yawThresh = deg2rad(getScalar(cfg, 'assocMaxYawRateDegPerSec', inf));
    pitchThresh = deg2rad(getScalar(cfg, 'assocMaxPitchRateDegPerSec', inf));

    for k = 1:nFrame
        meas(k).yawRate = yawRate(k);
        meas(k).pitchRate = pitchRate(k);
        meas(k).isYawing = yawRate(k) > yawThresh;
        meas(k).isPitching = pitchRate(k) > pitchThresh;
    end
end

function value = getScalar(s, name, defaultValue)
    value = defaultValue;
    if isfield(s, name) && ~isempty(s.(name))
        value = s.(name);
    end
end

function angle = wrapAngle(angle)
    angle = atan2(sin(angle), cos(angle));
end
