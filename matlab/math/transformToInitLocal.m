function [TLocalAll, timestamp] = transformToInitLocal(TGlobalAll, TinitGlobal, mode)
% transformToInitLocal Convert global SE(3) poses to an init-pose local frame.
%
% 입력:
%   TGlobalAll:
%     4 x 4 x N SE(3) matrices. Prefer converting raw row logs once with
%     pose12ToSE3 before calling this function.
%     Nx12/Nx13 row pose logs are also accepted for compatibility.
%
%   TinitGlobal:
%     4 x 4 SE(3) reference pose in global/UTM frame. If omitted, the first
%     pose in TGlobalAll is used.
%
%   mode:
%     'full' : full SE(3) first-pose 기준
%              Tlocal_k = inv(T0) * Tk
%
%     'yaw'  : road-plane/mapLocal 기준
%              XY는 Rz(yaw0)' * (xy - xy0)
%              Z는 z - z0
%              R은 Rz(yaw0)' * Rk
%
% 출력:
%   TLocalAll 4 x 4 x N local SE(3) matrices
%   timestamp Nx1 timestamp, or [] when input has no timestamp

    if nargin == 2 && (ischar(TinitGlobal) || isstring(TinitGlobal))
        mode = TinitGlobal;
        TinitGlobal = [];
    end

    if nargin < 3 || isempty(mode)
        mode = 'full';
    end
    mode = lower(char(mode));

    timestamp = [];
    if isnumeric(TGlobalAll) && size(TGlobalAll, 1) == 4 && size(TGlobalAll, 2) == 4
        TAll = TGlobalAll;
    else
        [TAll, timestamp] = pose12ToSE3(TGlobalAll);
    end

    N = size(TAll, 3);
    if nargin < 2 || isempty(TinitGlobal)
        TinitGlobal = TAll(:, :, 1);
    end
    if ~isequal(size(TinitGlobal), [4, 4])
        error('TinitGlobal must be a 4x4 SE(3) matrix.');
    end

    switch mode
        case 'full'
            Tref = TinitGlobal;

        case 'yaw'
            R0 = TinitGlobal(1:3, 1:3);
            t0 = TinitGlobal(1:3, 4);
            yaw0 = atan2(R0(2,1), R0(1,1));

            Rz0 = [cos(yaw0), -sin(yaw0), 0;
                   sin(yaw0),  cos(yaw0), 0;
                   0,          0,         1];

            Tref = eye(4);
            Tref(1:3, 1:3) = Rz0;
            Tref(1:3, 4) = t0;

        otherwise
            error("mode must be 'full' or 'yaw'.");
    end

    TrefInv = invSE3_local(Tref);
    TLocalAll = repmat(eye(4), [1, 1, N]);

    for k = 1:N
        TLocalAll(:, :, k) = TrefInv * TAll(:, :, k);
    end
end

function Tinv = invSE3_local(T)
    R = T(1:3, 1:3);
    t = T(1:3, 4);

    Tinv = eye(4);
    Tinv(1:3, 1:3) = R';
    Tinv(1:3, 4) = -R' * t;
end
