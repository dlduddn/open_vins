function xTrue = initialAlignmentTruth(cfg, mapRefPose, queryPose)
%INITIALALIGNMENTTRUTH True query pose in the map-initialization frame.
%
% buildMap() defines its local frame with mapRefPose plus
% cfg.mapInitXError/YError/YawError. This function expresses queryPose in
% that perturbed frame, so it can be compared with x0Align.

    if nargin < 3 || isempty(queryPose)
        queryPose = mapRefPose;
    end

    Rref0 = mapRefPose(1:3, 1:3);
    yawRef0 = atan2(Rref0(2, 1), Rref0(1, 1));
    tRef0 = mapRefPose(1:2, 4);

    Rquery = queryPose(1:3, 1:3);
    yawQuery = atan2(Rquery(2, 1), Rquery(1, 1));
    tQuery = queryPose(1:2, 4);

    dx = getScalar(cfg, 'mapInitXError', 0.0);
    dy = getScalar(cfg, 'mapInitYError', 0.0);
    dyaw = getScalar(cfg, 'mapInitYawError', 0.0);

    yawRef = yawRef0 + dyaw;
    Rref = [cos(yawRef), -sin(yawRef);
            sin(yawRef),  cos(yawRef)];
    tRef = tRef0 + [dx; dy];

    xTrue = [Rref' * (tQuery - tRef); wrapAngle(yawQuery - yawRef)];
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
