function statesOut = recenterStatesByPose(statesIn, xPose)
%RECENTERSTATESBYPOSE Express SE(2) states in the body frame defined by xPose.
%
% xPose maps the new frame into the current map frame:
%   p_map = R(yaw) * p_new + [x; y]
%
% statesIn is 3xN [x; y; yaw] in the current map frame.

    if isempty(statesIn)
        statesOut = statesIn;
        return;
    end
    if size(statesIn, 1) ~= 3
        error('statesIn must be a 3xN [x; y; yaw] matrix.');
    end

    xPose = xPose(:);
    if numel(xPose) ~= 3 || any(~isfinite(xPose))
        error('xPose must be a finite [x; y; yaw] vector.');
    end

    c = cos(xPose(3));
    s = sin(xPose(3));
    R = [c, -s; s, c];

    statesOut = statesIn;
    valid = all(isfinite(statesIn), 1);
    statesOut(1:2, valid) = R' * (statesIn(1:2, valid) - xPose(1:2));
    statesOut(3, valid) = wrapAngle(statesIn(3, valid) - xPose(3));
end

function angle = wrapAngle(angle)
    angle = atan2(sin(angle), cos(angle));
end
