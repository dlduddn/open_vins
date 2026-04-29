function [interpTimestamps, interpPoses] = genRef(cfg, queryTimestamps)
% interpolateGTPose  Query timestamp에 맞춰 GT pose를 SE(3) 보간
%
% 입력:
%   cfg.imageDir     - 이미지 폴더 경로 (fallback: 파일명 = 나노초 ROS timestamp)
%   cfg.gtDir        - GT pose CSV 경로
%   queryTimestamps  - optional. 보간할 timestamp. seconds 또는 nanoseconds 허용.
%
% 출력:
%   interpTimestamps - (N,1) 보간에 사용된 timestamp [ns]
%   interpPoses      - (N,12) 보간된 pose [P(0,0)~P(2,3)]

%% 1. Query timestamp 확정
queryTimestamps = sort(queryTimestamps(:));

%% 2. GT pose 로드
gtData = readmatrix(cfg.gtDir);
gtTimestamps = gtData(:, 1) * 1e-9; % ns -> sec
gtPoses      = gtData(:, 2:13);
fprintf('[GT] count = %d, range = %.3f ~ %.3f sec\n', ...
    size(gtData, 1), ...
    gtTimestamps(1) - queryTimestamps(1), ...
    gtTimestamps(end) - queryTimestamps(1));

%% 3. Query/GT 시간 범위 정렬
% GT 시작 이전 query 구간이 있으면, 첫 GT pose를 복제하여 앞쪽 보간용 GT를 생성
preGtQueryMask = queryTimestamps < gtTimestamps(1);
preGtQueryTimestamps = queryTimestamps(preGtQueryMask);

if ~isempty(preGtQueryTimestamps)
    numPrependedGt = numel(preGtQueryTimestamps);
    prependedGtPoses = repmat(gtPoses(1, :), numPrependedGt, 1);

    gtTimestamps = [preGtQueryTimestamps; gtTimestamps];
    gtPoses      = [prependedGtPoses; gtPoses];

    % fprintf('[GT] prepended first pose for %d pre-GT query timestamps\n', numPrependedGt);
end

% GT 종료 이후 query timestamp는 제거
postGtQueryMask = queryTimestamps > gtTimestamps(end);
numDroppedQueries = sum(postGtQueryMask);

if numDroppedQueries > 0
    queryTimestamps = queryTimestamps(~postGtQueryMask);
    % fprintf('[%s] removed %d post-GT queries, kept %d\n', ...
    %     queryLabel, numDroppedQueries, numel(queryTimestampsNsec));
end

%% 4. GT 무결성 점검 및 누락 pose 보간
% GT 기준 샘플링 주기와 결측 판정 임계값 설정
nominalGtDtNsec   = 0.01 * 1e9;             % 기준 샘플 간격(0.01초)을 ns 단위로 변환
gtGapThresholdNsec = 2 * nominalGtDtNsec;   % 인접 GT timestamp 간격이 0.02초를 초과하면 누락으로 판정

gtDtNsec   = diff(gtTimestamps);                         % 인접 GT timestamp 간 시간 차이 계산
gapStartIdx = find(gtDtNsec > gtGapThresholdNsec);      % 누락 구간 시작 인덱스 탐지

if ~isempty(gapStartIdx)
    % fprintf('[GT-Integrity] detected %d gap segments\n', numel(gapStartIdx));

    interpGapTimestamps = [];
    interpGapPoses      = [];

    % 각 누락 구간에 대해 translation은 선형 보간, rotation은 SLERP 수행
    for g = 1:numel(gapStartIdx)
        k = gapStartIdx(g);

        % 양 끝 GT timestamp 추출
        tStart = gtTimestamps(k);
        tEnd   = gtTimestamps(k + 1);

        % 누락 구간 내부에 삽입할 보간 timestamp 생성 (양 끝 제외)
        numFill = round((tEnd - tStart) / nominalGtDtNsec) - 1;
        if numFill < 1
            continue;
        end
        fillTimestamps = tStart + nominalGtDtNsec * (1:numFill)';

        % 양 끝 GT pose의 rotation/translation 추출
        [rotStart, transStart] = unpackPoseRow(gtPoses(k, :));
        [rotEnd,   transEnd]   = unpackPoseRow(gtPoses(k + 1, :));

        quatStart = quaternion(rotm2quat(rotStart));
        quatEnd   = quaternion(rotm2quat(rotEnd));

        % 누락 구간 내부에 삽입할 보간 pose 생성 (양 끝 제외)
        fillPoses = zeros(numFill, 12);
        for j = 1:numFill
            alpha = (fillTimestamps(j) - tStart) / (tEnd - tStart);

            % Translation 선형 보간
            transInterp = (1 - alpha) * transStart + alpha * transEnd;

            % Rotation SLERP
            rotInterp = rotmat(slerp(quatStart, quatEnd, alpha), 'point');

            fillPoses(j, :) = packPoseRow(rotInterp, transInterp);
        end

        interpGapTimestamps = [interpGapTimestamps; fillTimestamps]; %#ok<AGROW>
        interpGapPoses      = [interpGapPoses; fillPoses]; %#ok<AGROW>
    end

    % 기존 GT에 보간 데이터를 병합한 뒤 timestamp 기준으로 재정렬
    if ~isempty(interpGapTimestamps)
        gtTimestamps = [gtTimestamps; interpGapTimestamps];
        gtPoses      = [gtPoses; interpGapPoses];
        [gtTimestamps, sortIdx] = sort(gtTimestamps);
        gtPoses = gtPoses(sortIdx, :);

        % fprintf('[GT-Integrity] inserted %d poses, total GT=%d\n', ...
            % numel(interpGapTimestamps), numel(gtTimestamps));
    end
else
    % fprintf('[GT-Integrity] no gaps detected\n');
end

%% 5. GT pose를 quaternion과 translation으로 변환
numGtPoses = size(gtPoses, 1);
gtQuat  = zeros(numGtPoses, 4);
gtTrans = zeros(numGtPoses, 3);

for i = 1:numGtPoses
    [rotMat, transVec] = unpackPoseRow(gtPoses(i, :));
    gtQuat(i, :)  = rotm2quat(rotMat);   % [w x y z]
    gtTrans(i, :) = transVec;
end

%% 6. 보간에 사용할 query timestamp 확정 및 최종 pose 보간
interpTimestamps = queryTimestamps';
fprintf('[Interp] using %d query timestamps\n', numel(interpTimestamps));

assert(all(diff(gtTimestamps) > 0), 'gtTimestamps must be strictly increasing');

numInterp = numel(interpTimestamps);
interpPoses = zeros(numInterp, 12);

for i = 1:numInterp
    queryTimestamp = interpTimestamps(i);

    % GT 범위 밖 클램핑
    if queryTimestamp <= gtTimestamps(1)
        interpPoses(i, :) = gtPoses(1, :);
        continue;
    elseif queryTimestamp >= gtTimestamps(end)
        interpPoses(i, :) = gtPoses(end, :);
        continue;
    end

    leftIdx = find(gtTimestamps <= queryTimestamp, 1, 'last');

    t0 = gtTimestamps(leftIdx);
    t1 = gtTimestamps(leftIdx + 1);

    alpha = (queryTimestamp - t0) / (t1 - t0);
    alpha = max(0, min(1, alpha));

    % Translation 선형 보간
    transInterp = (1 - alpha) * gtTrans(leftIdx, :) + alpha * gtTrans(leftIdx + 1, :);

    % Rotation SLERP
    quat0 = quaternion(gtQuat(leftIdx, :));
    quat1 = quaternion(gtQuat(leftIdx + 1, :));
    quatInterp = slerp(quat0, quat1, alpha);
    rotInterp = rotmat(quatInterp, 'point');

    interpPoses(i, :) = packPoseRow(rotInterp, transInterp);
end

fprintf('[Interp] completed %d interpolated poses\n', numInterp);

end

function [rotMat, transVec] = unpackPoseRow(poseRow)
    % unpackPoseRow  1x12 pose row를 rotation matrix와 translation vector로 분해
    rotMat = [poseRow(1:3); poseRow(5:7); poseRow(9:11)];
    transVec = [poseRow(4), poseRow(8), poseRow(12)];
end

function poseRow = packPoseRow(rotMat, transVec)
    % packPoseRow  rotation matrix와 translation vector를 1x12 pose row로 결합
    poseRow = zeros(1, 12);
    poseRow(1:3)  = rotMat(1, :);
    poseRow(4)    = transVec(1);
    poseRow(5:7)  = rotMat(2, :);
    poseRow(8)    = transVec(2);
    poseRow(9:11) = rotMat(3, :);
    poseRow(12)   = transVec(3);
end
