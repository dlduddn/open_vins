function measurement = loadYaml(cfg, yamlName)
    % YAML의 bezier 섹션에서 [X Y width_per_point id]만 읽는다.
    %
    % 출력 measurement.points는 모든 bezier curve point를 하나로 합친 N x 4 행렬이다.
    % 각 row는 [pixel_x, pixel_y, width_per_point, curve_id] 형식이다.
    % measurement.curves에는 curve별 원본 pixel/width와 curve별 points를 따로 저장한다.
    measurement = struct();
    measurement.yamlName = char(yamlName);
    measurement.points = zeros(0, 4);
    measurement.curves = struct('id', {}, 'pixel', {}, 'width_per_point', {}, 'points', {});

    % yamlName이 없으면 매칭된 YAML이 없다는 뜻이므로 빈 measurement를 반환한다.
    if nargin < 2 || isempty(yamlName)
        return;
    end

    % yamlName은 확장자 없이 들어올 수도 있고, .yaml까지 포함해서 들어올 수도 있다.
    % 확장자가 없으면 기본적으로 .yaml을 붙인다.
    [~, name, ext] = fileparts(char(yamlName));
    if isempty(ext), ext = '.yaml'; end

    % cfg.yamlDir 아래에서 실제 YAML 파일 경로를 만든다.
    % 파일이 없으면 오류 대신 빈 measurement를 반환해 호출부가 자연스럽게 skip할 수 있게 한다.
    yamlPath = fullfile(cfg.yamlDir, [name ext]);
    if ~isfile(yamlPath)
        return;
    end

    % YAML 파일의 bezier curve 정보를 읽는다.
    measurement.yamlName = name;
    measurement.curves = parseBezierYaml(yamlPath);

    % curve별 pixel 좌표와 width_per_point를 [X Y width id] 형태로 결합한다.
    for id = 1:numel(measurement.curves)
        pixel = measurement.curves(id).pixel;
        widthPerPoint = measurement.curves(id).width_per_point(:);

        % 각 pixel point마다 width 값이 하나씩 있어야 한다.
        if numel(widthPerPoint) ~= size(pixel, 1)
            error('width_per_point length does not match point count in %s (curve id=%d)', yamlPath, id);
        end

        % curve id를 마지막 column에 붙여, 나중에 전체 point 행렬에서도 어떤 curve에서 왔는지 알 수 있게 한다.
        points = [pixel, widthPerPoint, repmat(id, size(pixel, 1), 1)];
        measurement.curves(id).id = id;
        measurement.curves(id).points = points;

        % 전체 curve의 point를 하나의 행렬로 누적한다.
        measurement.points = [measurement.points; points]; %#ok<AGROW>
    end
end

function curves = parseBezierYaml(yamlPath)
    % MATLAB 기본 환경에 YAML parser가 없을 수 있으므로 필요한 형식만 직접 파싱한다.
    % 이 함수는 bezier 섹션 안의 curve_type, pixel, width_per_point만 이해한다.

    % OS별 줄바꿈(CRLF/LF/CR)을 모두 처리하기 위해 regexp로 line split한다.
    lines = regexp(fileread(yamlPath), '\r\n|\n|\r', 'split')';
    curves = struct('id', {}, 'pixel', {}, 'width_per_point', {}, 'points', {});

    % inBezier는 현재 읽는 위치가 YAML의 bezier: 섹션 안인지 나타낸다.
    inBezier = false;

    % curve는 현재 파싱 중인 bezier curve 하나를 임시로 담는 구조체이다.
    curve = [];
    i = 1;

    while i <= numel(lines)
        line = lines{i};
        text = strtrim(line);

        % bezier: 섹션을 만나면 이후 line부터 curve 정보를 읽기 시작한다.
        if strcmp(text, 'bezier:')
            inBezier = true;
            i = i + 1;
            continue;
        end

        % bezier 섹션 안에서 indentation 없는 새 top-level key를 만나면
        % bezier 섹션이 끝난 것으로 보고 파싱을 종료한다.
        if inBezier && ~isempty(regexp(line, '^[A-Za-z_][A-Za-z0-9_]*:\s*$', 'once'))
            break;
        end

        % bezier 섹션에 들어가기 전 line들은 모두 무시한다.
        if ~inBezier
            i = i + 1;
            continue;
        end

        % 새 bezier curve 항목을 만나면 이전 curve를 저장하고 새 curve 구조체를 시작한다.
        if startsWith(text, '- curve_type: bezier')
            if ~isempty(curve) && ~isempty(curve.pixel)
                curves(end + 1) = curve; %#ok<AGROW>
            end

            % points는 loadYaml에서 [X Y width id]로 채우므로 여기서는 빈 값으로 둔다.
            curve = struct('id', [], 'pixel', zeros(0, 2), ...
                'width_per_point', [], 'points', zeros(0, 4));
            i = i + 1;
            continue;
        end

        % curve_type을 만나기 전의 bezier 섹션 내부 line은 무시한다.
        if isempty(curve)
            i = i + 1;
            continue;
        end

        % pixel: 아래의 숫자 block을 읽어 [x1; y1; x2; y2; ...] 형태의 vector로 받는다.
        % 이후 2개씩 묶어 N x 2 [X Y] 행렬로 바꾼다.
        if strcmp(text, 'pixel:')
            [values, i] = readNumberBlock(lines, i + 1);
            if mod(numel(values), 2) ~= 0
                error('pixel field must contain X/Y pairs in %s', yamlPath);
            end
            curve.pixel = reshape(values, 2, []).';
            continue;
        elseif strcmp(text, 'width_per_point:')
            % width_per_point: 아래의 숫자 block은 point 개수와 같은 길이여야 한다.
            [curve.width_per_point, i] = readNumberBlock(lines, i + 1);
            continue;
        end

        % 관심 없는 field는 건너뛴다.
        i = i + 1;
    end

    % 파일 끝에서 마지막 curve가 아직 curves에 들어가지 않았으면 저장한다.
    if ~isempty(curve) && ~isempty(curve.pixel)
        curves(end + 1) = curve;
    end
end

function [values, nextIdx] = readNumberBlock(lines, startIdx)
    % YAML list item 형태로 이어진 숫자 block을 읽는다.
    % 지원하는 형태:
    %   - 123.4
    %   - -123.4
    %   - 1.2e-3
    number = '[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?';
    values = [];
    nextIdx = startIdx;

    while nextIdx <= numel(lines)
        % 현재 line이 "- 숫자" 형식이면 숫자를 읽고, 아니면 block이 끝난 것으로 본다.
        % 첫 번째 "- "는 YAML list marker이고, 그 뒤의 token을 숫자로 읽는다.
        token = regexp(lines{nextIdx}, ['^\s*-\s*(?:-\s*)?(' number ')\s*$'], 'tokens', 'once');
        if isempty(token)
            break;
        end

        % 문자열 token을 double로 바꿔 누적한다.
        values = [values; str2double(token{1})]; %#ok<AGROW>
        nextIdx = nextIdx + 1;
    end
end
