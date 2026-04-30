function [intrinsics, resolution, info] = loadKalibrCameraConfig(yamlPath, camName)
%LOADKALIBRCAMERACONFIG Read intrinsics/resolution from a Kalibr camera YAML.

    if nargin < 2 || isempty(camName)
        camName = 'cam0';
    end
    camName = char(camName);

    intrinsics = [];
    resolution = [];
    info = struct('yamlPath', yamlPath, 'camName', camName);

    if nargin < 1 || isempty(yamlPath) || ~isfile(yamlPath)
        error('Kalibr camera YAML not found: %s', string(yamlPath));
    end

    text = fileread(yamlPath);
    lines = splitlines(string(text));

    startIdx = find(lines == string(camName) + ":", 1, 'first');
    if isempty(startIdx)
        error('Camera section "%s" not found in %s', camName, yamlPath);
    end

    endIdx = numel(lines);
    for k = startIdx + 1:numel(lines)
        if ~isempty(regexp(char(lines(k)), '^cam[0-9]+:\s*$', 'once'))
            endIdx = k - 1;
            break;
        end
    end

    section = strjoin(lines(startIdx:endIdx), newline);
    intrinsics = parseBracketVector(section, 'intrinsics', 4);
    resolution = parseBracketVector(section, 'resolution', 2);

    info.intrinsics = intrinsics;
    info.resolution = resolution;
end

function values = parseBracketVector(section, key, expectedCount)
    pattern = string(key) + "\s*:\s*\[([^\]]+)\]";
    token = regexp(section, char(pattern), 'tokens', 'once');
    if isempty(token)
        error('Key "%s" not found in selected Kalibr camera section.', key);
    end

    values = regexp(token{1}, '[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?', 'match');
    values = str2double(values);
    if numel(values) < expectedCount || any(~isfinite(values(1:expectedCount)))
        error('Key "%s" does not contain %d finite numeric values.', key, expectedCount);
    end
    values = values(1:expectedCount);
end
