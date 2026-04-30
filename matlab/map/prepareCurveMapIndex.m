function mapIndex = prepareCurveMapIndex(cfg, mapDB)
%PREPARECURVEMAPINDEX Cache valid map rows and nearest-neighbor structures.
%
% curveMapLikelihood is called many times during initialization and PF.
% Building the valid map arrays and KD-tree once avoids repeating that work
% for every candidate pose.

    if nargin < 2 || isempty(mapDB) || size(mapDB, 2) < 2
        mapIndex = emptyMapIndex();
        return;
    end

    valid = all(isfinite(mapDB(:, 1:2)), 2);
    mapXY = mapDB(valid, 1:2);
    if size(mapDB, 2) >= 4
        mapWidth = mapDB(valid, 4);
    else
        mapWidth = NaN(size(mapXY, 1), 1);
    end

    mapIndex = emptyMapIndex();
    mapIndex.xy = mapXY;
    mapIndex.width = mapWidth;
    mapIndex.hasKdTree = false;
    mapIndex.searcher = [];

    if isempty(mapXY)
        return;
    end

    mapIndex.maxRoadWidth = maxFinite(mapWidth, getScalar(cfg, 'assocDefaultRoadWidth', 6.0));

    useKdTree = getLogical(cfg, 'assocUseMapKdTree', true);
    if useKdTree && exist('createns', 'file') == 2 && exist('knnsearch', 'file') == 2
        mapIndex.searcher = createns(mapXY, 'NSMethod', 'kdtree');
        mapIndex.hasKdTree = true;
    end
end

function mapIndex = emptyMapIndex()
    mapIndex = struct( ...
        'xy', zeros(0, 2), ...
        'width', zeros(0, 1), ...
        'maxRoadWidth', NaN, ...
        'hasKdTree', false, ...
        'searcher', []);
end

function value = maxFinite(values, defaultValue)
    values = values(isfinite(values) & values > 0);
    if isempty(values)
        value = defaultValue;
    else
        value = max(values);
    end
end

function value = getScalar(s, name, defaultValue)
    value = defaultValue;
    if isfield(s, name) && ~isempty(s.(name))
        value = s.(name);
    end
end

function value = getLogical(s, name, defaultValue)
    value = defaultValue;
    if isfield(s, name) && ~isempty(s.(name))
        value = logical(s.(name));
    end
end
