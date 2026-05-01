function angle = wrapAngle(angle)
%WRAPANGLE Normalize angles to [-pi, pi].

    angle = atan2(sin(angle), cos(angle));
end
