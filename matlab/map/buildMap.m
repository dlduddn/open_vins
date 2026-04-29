function [map] = buildMap(cfg, firstPose)
%BUILDSDMAPPOINTCLOUD  SHP 도로중심선을 리샘플링하여 첫 프레임 좌표 2D 점군 생성
%
% SD Map (Shapefile 도로중심선)을 균일 리샘플링하고,
% firstPose 기준 좌표계로 정규화하여 2D 점군을 반환한다.
%
% 전체 처리 흐름
%   1) SHP 로드
%   2) 필요 시 EPSG:5179 -> UTM 변환
%   3) 도로유형 필터링
%   4) 각 중심선을 일정 간격으로 리샘플링
%   5) 첫 프레임 좌표계(yaw + XY) 기준으로 정규화
%
% 입력:
%   sdCfg     - 설정 구조체:
%     .shpPath       - Shapefile 경로
%     .shpCRS        - 'epsg5179' 또는 'utm52n'
%     .utmZone       - UTM zone (default: 52)
%     .ds            - 리샘플링 간격 m (default: 0.5)
%     .showRoadTypes - 도로유형 필터 (default: {'RDD000'})
%     .mapInitXError   - reference pose x 오차, global/UTM east [m]
%     .mapInitYError   - reference pose y 오차, global/UTM north [m]
%     .mapInitYawError - reference pose yaw 오차 [rad]
%   firstPose - 4x4 SE(3) reference pose in global/UTM frame
%
% 출력:
%   sdPoints - (P,4) 점군 + 속성 [x, y, RDLN, RVWD] (첫 프레임 좌표계)

% ---------- 기본 설정값 보완 ----------
if ~isfield(cfg, 'utmZone'),       cfg.utmZone = 52;          end
if ~isfield(cfg, 'ds'),            cfg.ds = 0.5;              end
if ~isfield(cfg, 'showRoadTypes'), cfg.showRoadTypes = {'RDD000'}; end
if ~isfield(cfg, 'mapInitXError'),   cfg.mapInitXError = 0.0;   end
if ~isfield(cfg, 'mapInitYError'),   cfg.mapInitYError = 0.0;   end
if ~isfield(cfg, 'mapInitYawError'), cfg.mapInitYawError = 0.0; end

%% SHP 로드 및 좌표 변환
fprintf('SD Map: Shapefile 로딩: %s\n', cfg.shpPath);
S = shaperead(cfg.shpPath, 'UseGeoCoords', false);

% 입력 SHP가 EPSG:5179라면, 이후 연산 편의를 위해 UTM으로 변환
if strcmp(cfg.shpCRS, 'epsg5179')
    fprintf('  EPSG:5179 → UTM Zone %dN 변환 중...\n', cfg.utmZone);
    for i = 1:numel(S)
        % shapefile 좌표는 NaN으로 선분이 끊기는 경우가 있으므로
        % 유효 좌표만 뽑아서 변환하고, 다시 원래 크기로 복원한다.
        x5179 = S(i).X(:);
        y5179 = S(i).Y(:);
        valid = ~isnan(x5179) & ~isnan(y5179);

        [lat, lon] = korea2000ToLatlon(x5179(valid), y5179(valid));
        [utmE, utmN] = latlon2utm(lat, lon, cfg.utmZone);

        xOut = nan(size(x5179));  yOut = nan(size(y5179));
        xOut(valid) = utmE;       yOut(valid) = utmN;
        S(i).X = xOut';           S(i).Y = yOut';
    end
    fprintf('  변환 완료.\n');
end

% 도로유형 필터
if ~isempty(cfg.showRoadTypes) && isfield(S, 'RDDV')
    rddvAll = {S.RDDV};
    S = S(ismember(rddvAll, cfg.showRoadTypes));
    fprintf('  도로유형 필터 [%s]: %d segments\n', ...
        strjoin(cfg.showRoadTypes, ', '), numel(S));
end

%% 3. 중심선 균일 리샘플링 + 점별 속성 부여
fprintf('  리샘플링 (ds=%.2f m)...\n', cfg.ds);

allPts = [];
allLaneCounts = [];
allRoadWidths = [];

for i = 1:numel(S)
    x = S(i).X(:);
    y = S(i).Y(:);

    % NaN separator 제거
    valid = ~isnan(x) & ~isnan(y);
    xv = x(valid);
    yv = y(valid);

    % 점이 2개 미만이면 선분 길이를 정의할 수 없으므로 건너뜀
    if numel(xv) < 2
        continue;
    end

    % 선분 길이와 누적 거리 계산
    segLen  = sqrt(diff(xv).^2 + diff(yv).^2);
    cumDist = [0; cumsum(segLen)];
    totalLen = cumDist(end);

    if totalLen < cfg.ds
        % 너무 짧은 선분은 원본 점 그대로 사용
        pts = [xv, yv];
    else
        % 0 ~ totalLen 구간을 ds 간격으로 샘플링
        sQuery = (0:cfg.ds:totalLen)';

        % 마지막 점이 빠지지 않도록 totalLen을 한 번 더 포함
        if sQuery(end) < totalLen
            sQuery = [sQuery; totalLen]; %#ok<AGROW>
        end

        % 누적 거리 축 기준 선형 보간
        pts = [interp1(cumDist, xv, sQuery, 'linear'), ...
               interp1(cumDist, yv, sQuery, 'linear')];
    end

    % 세그먼트 속성 읽기
    laneCount = getNumericField(S, i, 'RDLN');
    roadWidth = getNumericField(S, i, 'RVWD');

    % 샘플링된 각 점에 세그먼트 속성 반복 부여
    nPts = size(pts, 1);

    allPts = [allPts; pts]; %#ok<AGROW>
    allLaneCounts = [allLaneCounts; repmat(laneCount, nPts, 1)]; %#ok<AGROW>
    allRoadWidths = [allRoadWidths; repmat(roadWidth, nPts, 1)]; %#ok<AGROW>
end

fprintf('  UTM 점군: %d points\n', size(allPts, 1));

%% 4. 첫 프레임 로컬 좌표계로 정규화
% 이 문제는 2D road-plane 정합이므로 z/roll/pitch는 무시하고
% reference pose의 yaw와 XY 위치만 사용해 local map으로 변환한다.
T0 = firstPose;
R0 = T0(1:3, 1:3);
t0 = T0(1:2, 4);
t0 = t0 + [cfg.mapInitXError; cfg.mapInitYError];

% 첫 프레임 yaw 추출
yaw0 = atan2(R0(2,1), R0(1,1));
yaw0 = yaw0 + cfg.mapInitYawError;
R0_2d = [cos(yaw0), -sin(yaw0);
         sin(yaw0),  cos(yaw0)];

% UTM(E, N) -> 첫 프레임 2D 로컬 좌표계
sdXY = (R0_2d' * (allPts' - t0))';

% 점좌표 + 속성 직접 결합
map = [sdXY, allLaneCounts, allRoadWidths];

fprintf('  SD Map 점군 완료: %d points (첫 프레임 좌표계)\n', size(map, 1));
fprintf('  범위 X: [%.1f, %.1f] m,  Y: [%.1f, %.1f] m\n', ...
    min(map(:,1)), max(map(:,1)), ...
    min(map(:,2)), max(map(:,2)));

% 속성 요약 출력
validLane = ~isnan(map(:,3));
validWidth = ~isnan(map(:,4));
fprintf('  RDLN 부여 점 수: %d / %d\n', sum(validLane), size(map,1));
fprintf('  RVWD 부여 점 수: %d / %d\n', sum(validWidth), size(map,1));

end

%% ====================== LOCAL FUNCTIONS ======================
function [E, N] = latlon2utm(lat, lon, zone)
    a = 6378137.0;
    f = 1 / 298.257223563;
    e2 = 2*f - f^2;
    e_prime2 = e2 / (1 - e2);
    k0 = 0.9996;
    lon0 = (zone - 1) * 6 - 180 + 3;

    latR = deg2rad(lat);  lonR = deg2rad(lon);  lon0R = deg2rad(lon0);
    sinLat = sin(latR);  cosLat = cos(latR);  tanLat = tan(latR);
    Nrad = a ./ sqrt(1 - e2 * sinLat.^2);
    T = tanLat.^2;
    C = e_prime2 * cosLat.^2;
    A = cosLat .* (lonR - lon0R);

    M = a * ((1 - e2/4 - 3*e2^2/64 - 5*e2^3/256) .* latR ...
           - (3*e2/8 + 3*e2^2/32 + 45*e2^3/1024) .* sin(2*latR) ...
           + (15*e2^2/256 + 45*e2^3/1024) .* sin(4*latR) ...
           - (35*e2^3/3072) .* sin(6*latR));

    E = k0 * Nrad .* (A + (1-T+C).*A.^3/6 ...
        + (5-18*T+T.^2+72*C-58*e_prime2).*A.^5/120) + 500000;
    N = k0 * (M + Nrad .* tanLat .* (A.^2/2 + (5-T+9*C+4*C.^2).*A.^4/24 ...
        + (61-58*T+T.^2+600*C-330*e_prime2).*A.^6/720));
    N(lat < 0) = N(lat < 0) + 10000000;
end

function [lat, lon] = korea2000ToLatlon(easting, northing)
    a = 6378137.0;
    f = 1 / 298.257222101;
    k0 = 0.9996;
    E0 = 1000000;  N0 = 2000000;
    lat0 = 38 * pi / 180;
    lon0 = 127.5 * pi / 180;

    e2 = 2*f - f^2;
    e1 = (1 - sqrt(1-e2)) / (1 + sqrt(1-e2));

    x = easting - E0;
    y = northing - N0;

    M0 = meridianArc(lat0, a, e2);
    M  = M0 + y / k0;
    mu = M / (a * (1 - e2/4 - 3*e2^2/64 - 5*e2^3/256));

    lat1 = mu + (3*e1/2 - 27*e1^3/32)*sin(2*mu) ...
              + (21*e1^2/16 - 55*e1^4/32)*sin(4*mu) ...
              + (151*e1^3/96)*sin(6*mu) ...
              + (1097*e1^4/512)*sin(8*mu);

    N1 = a ./ sqrt(1 - e2 * sin(lat1).^2);
    T1 = tan(lat1).^2;
    C1 = (e2 / (1-e2)) * cos(lat1).^2;
    R1 = a * (1 - e2) ./ (1 - e2 * sin(lat1).^2).^1.5;
    D  = x ./ (N1 * k0);

    lat = lat1 - (N1 .* tan(lat1) ./ R1) .* ...
          (D.^2/2 - (5 + 3*T1 + 10*C1 - 4*C1.^2 - 9*e2/(1-e2)) .* D.^4/24 ...
          + (61 + 90*T1 + 298*C1 + 45*T1.^2 - 252*e2/(1-e2) - 3*C1.^2) .* D.^6/720);

    lon = lon0 + (D - (1 + 2*T1 + C1) .* D.^3/6 ...
          + (5 - 2*C1 + 28*T1 - 3*C1.^2 + 8*e2/(1-e2) + 24*T1.^2) .* D.^5/120) ./ cos(lat1);

    lat = lat * 180 / pi;
    lon = lon * 180 / pi;
end

function M = meridianArc(phi, a, e2)
    M = a * ((1 - e2/4 - 3*e2^2/64 - 5*e2^3/256) * phi ...
           - (3*e2/8 + 3*e2^2/32 + 45*e2^3/1024) * sin(2*phi) ...
           + (15*e2^2/256 + 45*e2^3/1024) * sin(4*phi) ...
           - (35*e2^3/3072) * sin(6*phi));
end

function value = getNumericField(S, idx, fieldName)
    % getNumericField
    % struct array의 field를 숫자로 읽어온다.
    % 필드가 없거나 숫자 변환이 불가능하면 NaN 반환.
    value = NaN;
    
    if ~isfield(S, fieldName)
        return;
    end
    
    raw = S(idx).(fieldName);
    
    if isnumeric(raw)
        if isempty(raw)
            value = NaN;
        else
            value = double(raw(1));
        end
    elseif ischar(raw) || isstring(raw)
        value = str2double(raw);
    else
        value = NaN;
    end
end
