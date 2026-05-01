function visualizePfMapMatching(cfg, mapIndex, t, len, xi, wi, xhat, neff, ...
        used, resampled, usableCount, meanLogL, matchViz)
%VISUALIZEPFMAPMATCHING Live diagnostic view for PF curve-map matching.

    if ~getLogical(cfg, 'pfVisualizeMapMatching', false)
        return;
    end

    interval = max(1, round(getScalar(cfg, 'pfVisualizeInterval', 1)));
    if ~(t == 1 || t == len || mod(t, interval) == 0)
        return;
    end

    persistent fig
    if isempty(fig) || ~ishandle(fig)
        fig = figure('Name', 'PF Map Matching', 'NumberTitle', 'off', 'Color', 'w');
    else
        figure(fig);
        clf(fig);
    end

    if isempty(fig) || ~ishandle(fig)
        return;
    end

    mapXY = getMapXY(mapIndex);
    score = getScoreVector(matchViz, wi);
    bestIdx = selectBestParticle(wi, score);
    center = chooseViewCenter(xhat, xi, bestIdx);
    viewHalfWidth = getScalar(cfg, 'pfVisualizeMapWindow', 45.0);

    axMap = subplot(2, 2, [1, 3], 'Parent', fig);
    hold(axMap, 'on');
    plotMapWindow(axMap, mapXY, center, viewHalfWidth);
    plotParticles(axMap, xi, wi, score, cfg);
    plotFrameCurves(axMap, matchViz, xhat, [0.0, 0.25, 0.95], '-', 'Estimate curves');
    if bestIdx > 0
        plotFrameCurves(axMap, matchViz, xi(:, bestIdx), [0.90, 0.10, 0.05], '--', 'Best particle curves');
        drawPoseArrow(axMap, xi(:, bestIdx), [0.90, 0.10, 0.05], 'Best particle');
    end
    drawPoseArrow(axMap, xhat, [0.0, 0.25, 0.95], 'Weighted estimate');
    formatMapAxis(axMap, center, viewHalfWidth);
    title(axMap, makeMapTitle(t, len, used, resampled, usableCount, neff, meanLogL, bestIdx, score));
    legend(axMap, 'Location', 'best');

    axBody = subplot(2, 2, 2, 'Parent', fig);
    plotBodyCurves(axBody, cfg, matchViz);
    title(axBody, 'Query Bezier curves in body frame');

    axScore = subplot(2, 2, 4, 'Parent', fig);
    plotScorePanel(axScore, score, wi, matchViz, bestIdx);

    drawnow limitrate;
end

function mapXY = getMapXY(mapIndex)
    if isstruct(mapIndex) && isfield(mapIndex, 'xy')
        mapXY = mapIndex.xy;
    else
        mapXY = zeros(0, 2);
    end
    if size(mapXY, 2) >= 2
        mapXY = mapXY(:, 1:2);
        mapXY = mapXY(all(isfinite(mapXY), 2), :);
    else
        mapXY = zeros(0, 2);
    end
end

function score = getScoreVector(matchViz, wi)
    n = numel(wi);
    score = [];
    if isstruct(matchViz) && isfield(matchViz, 'safeLogL') && numel(matchViz.safeLogL) == n
        score = matchViz.safeLogL(:).';
    elseif isstruct(matchViz) && isfield(matchViz, 'logL') && numel(matchViz.logL) == n
        score = matchViz.logL(:).';
    end
end

function bestIdx = selectBestParticle(wi, score)
    bestIdx = 0;
    if ~isempty(score) && any(isfinite(score))
        finiteScore = score;
        finiteScore(~isfinite(finiteScore)) = -inf;
        [~, bestIdx] = max(finiteScore);
    elseif ~isempty(wi) && any(isfinite(wi))
        finiteWeight = wi;
        finiteWeight(~isfinite(finiteWeight)) = -inf;
        [~, bestIdx] = max(finiteWeight);
    end
end

function center = chooseViewCenter(xhat, xi, bestIdx)
    if numel(xhat) >= 2 && all(isfinite(xhat(1:2)))
        center = xhat(1:2).';
    elseif bestIdx > 0 && size(xi, 2) >= bestIdx && all(isfinite(xi(1:2, bestIdx)))
        center = xi(1:2, bestIdx).';
    elseif ~isempty(xi)
        finite = all(isfinite(xi(1:2, :)), 1);
        if any(finite)
            center = mean(xi(1:2, finite), 2).';
        else
            center = [0, 0];
        end
    else
        center = [0, 0];
    end
end

function plotMapWindow(ax, mapXY, center, viewHalfWidth)
    if isempty(mapXY)
        return;
    end

    if isfinite(viewHalfWidth) && viewHalfWidth > 0
        keep = abs(mapXY(:, 1) - center(1)) <= viewHalfWidth & ...
               abs(mapXY(:, 2) - center(2)) <= viewHalfWidth;
        mapView = mapXY(keep, :);
    else
        mapView = mapXY;
    end

    if isempty(mapView)
        mapView = mapXY;
    end

    scatter(ax, mapView(:, 1), mapView(:, 2), 4, [0.72, 0.72, 0.72], ...
        'filled', 'DisplayName', 'Map centerline');
end

function plotParticles(ax, xi, wi, score, cfg)
    if isempty(xi) || size(xi, 1) < 2
        return;
    end

    nParticle = size(xi, 2);
    maxParticles = max(1, round(getScalar(cfg, 'pfVisualizeMaxParticles', nParticle)));
    idx = 1:nParticle;
    if nParticle > maxParticles
        [~, order] = sort(wi(:), 'descend');
        idx = order(1:maxParticles);
    end

    weights = wi(idx);
    markerSize = 10 + 60 * normalizePositive(weights);
    if ~isempty(score) && numel(score) == nParticle
        colorValue = score(idx);
        colorValue(~isfinite(colorValue)) = minFinite(score, -100.0) - 10.0;
        colorLabel = 'log L';
    else
        colorValue = weights;
        colorLabel = 'weight';
    end

    scatter(ax, xi(1, idx), xi(2, idx), markerSize, colorValue, 'filled', ...
        'MarkerEdgeColor', 'none', 'DisplayName', 'Particles');
    colormap(ax, parula);
    cb = colorbar(ax);
    ylabel(cb, colorLabel);
end

function out = normalizePositive(values)
    values = values(:).';
    values(~isfinite(values) | values < 0) = 0;
    maxValue = max(values);
    if isempty(maxValue) || maxValue <= 0
        out = zeros(size(values));
    else
        out = values / maxValue;
    end
end

function plotFrameCurves(ax, matchViz, x, color, lineStyle, displayName)
    curves = getUsableCurves(matchViz);
    if isempty(curves) || numel(x) < 3 || any(~isfinite(x(1:3)))
        return;
    end

    firstVisible = true;
    for i = 1:numel(curves)
        xy = transformBodyToMap(curves(i).sampleXY, x);
        if firstVisible
            plot(ax, xy(:, 1), xy(:, 2), lineStyle, 'Color', color, ...
                'LineWidth', 1.7, 'DisplayName', displayName);
            firstVisible = false;
        else
            plot(ax, xy(:, 1), xy(:, 2), lineStyle, 'Color', color, ...
                'LineWidth', 1.7, 'HandleVisibility', 'off');
        end
    end
end

function curves = getUsableCurves(matchViz)
    curves = [];
    if ~isstruct(matchViz) || ~isfield(matchViz, 'frame') || isempty(matchViz.frame)
        return;
    end

    frame = matchViz.frame;
    if ~isfield(frame, 'curves') || isempty(frame.curves)
        return;
    end

    if isfield(frame, 'usableIdx') && ~isempty(frame.usableIdx)
        curves = frame.curves(frame.usableIdx);
    else
        curves = frame.curves;
    end
end

function outXY = transformBodyToMap(bodyXY, x)
    c = cos(x(3));
    s = sin(x(3));
    outXY = bodyXY * [c, s; -s, c] + x(1:2).';
end

function drawPoseArrow(ax, x, color, displayName)
    if numel(x) < 3 || any(~isfinite(x(1:3)))
        return;
    end

    arrowLen = 4.0;
    plot(ax, x(1), x(2), 'o', 'Color', color, 'MarkerFaceColor', color, ...
        'MarkerSize', 5, 'HandleVisibility', 'off');
    quiver(ax, x(1), x(2), arrowLen * cos(x(3)), arrowLen * sin(x(3)), 0, ...
        'Color', color, 'LineWidth', 1.7, 'MaxHeadSize', 1.2, ...
        'DisplayName', displayName);
end

function formatMapAxis(ax, center, viewHalfWidth)
    axis(ax, 'equal');
    grid(ax, 'on');
    xlabel(ax, 'x map-local [m]');
    ylabel(ax, 'y map-local [m]');

    if isfinite(viewHalfWidth) && viewHalfWidth > 0
        xlim(ax, center(1) + [-viewHalfWidth, viewHalfWidth]);
        ylim(ax, center(2) + [-viewHalfWidth, viewHalfWidth]);
    end
end

function titleText = makeMapTitle(t, len, used, resampled, usableCount, neff, meanLogL, bestIdx, score)
    if ~isempty(score) && bestIdx > 0 && numel(score) >= bestIdx
        bestLogL = score(bestIdx);
    else
        bestLogL = NaN;
    end

    titleText = sprintf(['PF map matching %d/%d | used=%d resample=%d | ' ...
        'curves=%d N_eff=%.1f meanLogL=%.2f bestLogL=%.2f'], ...
        t, len, used, resampled, usableCount, neff, meanLogL, bestLogL);
end

function plotBodyCurves(ax, cfg, matchViz)
    cla(ax);
    hold(ax, 'on');
    drawCameraFov(ax, cfg);

    curves = getUsableCurves(matchViz);
    for i = 1:numel(curves)
        xy = curves(i).sampleXY;
        plot(ax, xy(:, 1), xy(:, 2), '-', 'Color', [0.0, 0.25, 0.95], ...
            'LineWidth', 1.5);
        plot(ax, xy(1, 1), xy(1, 2), '.', 'Color', [0.0, 0.25, 0.95], ...
            'MarkerSize', 12);
    end

    plot(ax, 0, 0, 'ks', 'MarkerFaceColor', 'k', 'MarkerSize', 5);
    axis(ax, 'equal');
    grid(ax, 'on');
    xlabel(ax, 'x body-forward [m]');
    ylabel(ax, 'y body-left [m]');

    maxRange = getScalar(cfg, 'cameraFovMaxRange', getScalar(cfg, 'bevH', 120) * getScalar(cfg, 'resolution', 0.5));
    if isfinite(maxRange) && maxRange > 0
        xlim(ax, [-5, maxRange + 5]);
        ylim(ax, [-maxRange * 0.7, maxRange * 0.7]);
    end
end

function drawCameraFov(ax, cfg)
    gate = computeCameraFovGate(cfg);
    if ~isfield(gate, 'enabled') || ~gate.enabled || ~isfinite(gate.maxRange)
        return;
    end

    angle = linspace(-gate.rightAngle - gate.margin, gate.leftAngle + gate.margin, 80) + gate.yawOffset;
    outer = gate.maxRange * [cos(angle(:)), sin(angle(:))];
    innerRange = max(gate.minRange, 0.0);
    inner = innerRange * [cos(flipud(angle(:))), sin(flipud(angle(:)))];
    poly = [outer; inner; outer(1, :)];

    patch(ax, poly(:, 1), poly(:, 2), [0.95, 0.90, 0.70], ...
        'FaceAlpha', 0.25, 'EdgeColor', [0.75, 0.55, 0.15], ...
        'LineStyle', '-', 'DisplayName', 'Camera FOV');
end

function plotScorePanel(ax, score, wi, matchViz, bestIdx)
    cla(ax);
    hold(ax, 'on');

    if ~isempty(score)
        finiteScore = score;
        finiteScore(~isfinite(finiteScore)) = minFinite(score, -100.0) - 10.0;
        sortedScore = sort(finiteScore, 'descend');
        plot(ax, sortedScore, '.-', 'Color', [0.10, 0.25, 0.65], 'LineWidth', 1.0);
        ylabel(ax, 'log likelihood');
        title(ax, makeScoreTitle(matchViz, bestIdx));
    else
        sortedWeight = sort(wi(:), 'descend');
        plot(ax, sortedWeight, '.-', 'Color', [0.10, 0.25, 0.65], 'LineWidth', 1.0);
        ylabel(ax, 'weight');
        title(ax, 'Particle weights');
    end

    grid(ax, 'on');
    xlabel(ax, 'particle rank');
end

function titleText = makeScoreTitle(matchViz, bestIdx)
    matchedCurves = NaN;
    meanResidual = NaN;
    if isstruct(matchViz) && isfield(matchViz, 'likelihoodInfo') && ...
            isstruct(matchViz.likelihoodInfo) && bestIdx > 0
        info = matchViz.likelihoodInfo;
        if isfield(info, 'matchedCurveCount') && numel(info.matchedCurveCount) >= bestIdx
            matchedCurves = info.matchedCurveCount(bestIdx);
        end
        if isfield(info, 'meanResidual') && numel(info.meanResidual) >= bestIdx
            meanResidual = info.meanResidual(bestIdx);
        end
    end

    titleText = sprintf('Likelihood rank | best matched=%g residual=%.2f m', ...
        matchedCurves, meanResidual);
end

function value = minFinite(values, defaultValue)
    values = values(isfinite(values));
    if isempty(values)
        value = defaultValue;
    else
        value = min(values);
    end
end
