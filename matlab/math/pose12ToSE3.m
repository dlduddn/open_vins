function [TAll, timestamp] = pose12ToSE3(poseIn)
% pose12ToSE3 Convert row pose logs to 4x4xN SE(3) matrices.
%
% 입력:
%   poseIn:
%     Nx12 [r00 r01 r02 tx r10 r11 r12 ty r20 r21 r22 tz]
%     Nx13 [timestamp r00 r01 r02 tx r10 r11 r12 ty r20 r21 r22 tz]
%
% 출력:
%   TAll      4 x 4 x N SE(3) matrices
%   timestamp Nx1 timestamp, or [] when poseIn is Nx12

    if size(poseIn, 2) == 13
        timestamp = poseIn(:, 1);
        pose12 = poseIn(:, 2:13);
    elseif size(poseIn, 2) == 12
        timestamp = [];
        pose12 = poseIn;
    else
        error('poseIn must be Nx12 or Nx13.');
    end

    N = size(pose12, 1);
    TAll = repmat(eye(4), [1, 1, N]);

    TAll(1, 1:3, :) = reshape(pose12(:, 1:3).', 1, 3, N);
    TAll(1, 4, :)   = reshape(pose12(:, 4).', 1, 1, N);
    TAll(2, 1:3, :) = reshape(pose12(:, 5:7).', 1, 3, N);
    TAll(2, 4, :)   = reshape(pose12(:, 8).', 1, 1, N);
    TAll(3, 1:3, :) = reshape(pose12(:, 9:11).', 1, 3, N);
    TAll(3, 4, :)   = reshape(pose12(:, 12).', 1, 1, N);
end
