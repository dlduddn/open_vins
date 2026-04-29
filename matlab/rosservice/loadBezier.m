function meas = loadBezier(cfg, epoch)
% loadBezierMeasurements Load timestamp-aligned Bezier YAML measurements.
%
% 입력:
%   cfg   - loadConfig()에서 만든 설정 구조체. yamlDir, bevW, bevH,
%           resolution, u0, v0, curveType 필드를 사용한다.
%   epoch - 1 x N 또는 N x 1 camera timestamp [sec]
%
% 출력:
%   meas  - 1 x N struct array. 각 원소는 epoch(k)에 대응하는 YAML 측정치이며,
%           points는 pixel 좌표, pointsBody는 body 좌표 [x_forward y_left width id]이다.
%           YAML 파일이 존재하지만 bezier: []이면 hasYaml=true, hasPoints=false로 저장한다.

    epoch = epoch(:).';
    nFrame = numel(epoch);

    emptyMeas = struct( ...
        't_cam', [], ...
        'yamlName', '', ...
        'points', zeros(0, 4), ...
        'pointsBody', zeros(0, 4));
    meas = repmat(emptyMeas, 1, nFrame);

    nYaml = 0;
    for k = 1:nFrame
        yamlName = searchYaml(cfg, epoch(k));
        measurement = loadYaml(cfg, yamlName);
        hasYaml = ~isempty(measurement.yamlName);

        meas(k).t_cam = epoch(k);
        meas(k).yamlName = measurement.yamlName; 
        meas(k).points = measurement.points;
        meas(k).pointsBody = pixelPointsToBody(cfg, measurement.points);

        nYaml = nYaml + double(hasYaml);
    end

    fprintf('Bezier YAML loaded: %d / %d frames\n', nYaml, nFrame);

end

function pointsBody = pixelPointsToBody(cfg, pointsPixel)
    if isempty(pointsPixel)
        pointsBody = zeros(0, 4);
        return;
    end

    uv = pointsPixel(:, 1:2);
    widthPixel = pointsPixel(:, 3);
    curveId = pointsPixel(:, 4);

    % BEV pixel frame:
    %   u increases to image-right, v increases downward.
    % Body/map SE(2) frame used by buildMap/sir:
    %   x is forward, y is left.
    xBody = (cfg.v0 - uv(:, 2)) * cfg.resolution;
    yBody = (cfg.u0 - uv(:, 1)) * cfg.resolution;
    widthBody = widthPixel * cfg.resolution;

    pointsBody = [xBody, yBody, widthBody, curveId];
end
