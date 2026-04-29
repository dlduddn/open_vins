function visualizeInitialAlignment(mapDB, alignInfo, xEst, xTrue)
%VISUALIZEINITIALALIGNMENT Plot initial map alignment result.

    if nargin < 4
        xTrue = [];
    end
    xEst = xEst(:);
    xTrue = xTrue(:);

    if ~isfield(alignInfo, 'sourceXY') || isempty(alignInfo.sourceXY)
        warning('visualizeInitialAlignment:MissingSource', ...
            'alignInfo has no sourceXY. Skipping visualization.');
        return;
    end

    sourceXY = alignInfo.sourceXY;
    mapXY = mapDB(:, 1:2);
    mapXY = mapXY(all(isfinite(mapXY), 2), :);
    if isfield(alignInfo, 'mapCrop') && ~isempty(alignInfo.mapCrop)
        mapView = alignInfo.mapCrop;
    else
        mapView = mapXY;
    end

    priorXY = getAlignedPoints(alignInfo, 'priorAlignedXY', sourceXY, ...
        getField(alignInfo, 'prior', [0; 0; 0]));
    coarseXY = getAlignedPoints(alignInfo, 'coarseAlignedXY', sourceXY, ...
        getField(alignInfo, 'coarseX0', [0; 0; 0]));
    estXY = getAlignedPoints(alignInfo, 'alignedXY', sourceXY, xEst);

    hasTrue = numel(xTrue) == 3 && all(isfinite(xTrue(:)));
    if hasTrue
        trueXY = transformPoints(sourceXY, xTrue);
        posErr = norm(xEst(1:2) - xTrue(1:2));
        yawErrDeg = rad2deg(wrapAngle(xEst(3) - xTrue(3)));
    else
        trueXY = [];
        posErr = NaN;
        yawErrDeg = NaN;
    end

    figure('Name', 'Initial Map Alignment');

    subplot(2, 2, [1, 3]);
    plotAlignmentScene(mapView, priorXY, coarseXY, estXY, trueXY, xEst, xTrue, hasTrue);
    title(sprintf('Initial alignment frame %d', getField(alignInfo, 'frameIdx', 0)));

    subplot(2, 2, 2);
    zoomXY = [estXY; trueXY; coarseXY];
    plotAlignmentScene(mapView, priorXY, coarseXY, estXY, trueXY, xEst, xTrue, hasTrue);
    setZoomAround(zoomXY, 15.0);
    title('Zoomed alignment');

    subplot(2, 2, 4);
    plotPoseComparison(xEst, xTrue, hasTrue);
    if hasTrue
        title(sprintf('Error: %.2fm, %.2fdeg', posErr, yawErrDeg));
    else
        title('Estimated pose');
    end
end

function plotAlignmentScene(mapXY, priorXY, coarseXY, estXY, trueXY, xEst, xTrue, hasTrue)
    scatter(mapXY(:, 1), mapXY(:, 2), 3, [0.75, 0.75, 0.75], ...
        'filled', 'DisplayName', 'Map centerline');
    hold on;
    scatter(priorXY(:, 1), priorXY(:, 2), 14, [0.95, 0.55, 0.20], ...
        'filled', 'DisplayName', 'Prior');
    scatter(coarseXY(:, 1), coarseXY(:, 2), 14, [0.15, 0.65, 0.75], ...
        'filled', 'DisplayName', 'Coarse');
    scatter(estXY(:, 1), estXY(:, 2), 18, [0.05, 0.25, 0.90], ...
        'filled', 'DisplayName', 'Estimate');

    if hasTrue
        scatter(trueXY(:, 1), trueXY(:, 2), 18, [0.05, 0.55, 0.20], ...
            'filled', 'DisplayName', 'True');
    end

    drawPoseArrow(xEst, [0.05, 0.25, 0.90], 'Estimate pose');
    if hasTrue
        drawPoseArrow(xTrue, [0.05, 0.55, 0.20], 'True pose');
    end

    axis equal;
    grid on;
    xlabel('x [m]');
    ylabel('y [m]');
    legend('Location', 'best');
end

function plotPoseComparison(xEst, xTrue, hasTrue)
    estVals = [xEst(1), xEst(2), rad2deg(xEst(3))];
    if hasTrue
        trueVals = [xTrue(1), xTrue(2), rad2deg(xTrue(3))];
        bar([trueVals; estVals]');
        legend('True', 'Estimate', 'Location', 'best');
    else
        bar(estVals');
        legend('Estimate', 'Location', 'best');
    end

    grid on;
    set(gca, 'XTickLabel', {'x [m]', 'y [m]', 'yaw [deg]'});
end

function drawPoseArrow(x, color, name)
    len = 3.0;
    quiver(x(1), x(2), len * cos(x(3)), len * sin(x(3)), 0, ...
        'Color', color, 'LineWidth', 2.0, 'MaxHeadSize', 1.0, ...
        'DisplayName', name);
end

function setZoomAround(xy, margin)
    xy = xy(all(isfinite(xy), 2), :);
    if isempty(xy)
        return;
    end

    xlim([min(xy(:, 1)) - margin, max(xy(:, 1)) + margin]);
    ylim([min(xy(:, 2)) - margin, max(xy(:, 2)) + margin]);
end

function xy = getAlignedPoints(info, fieldName, sourceXY, x)
    if isfield(info, fieldName) && ~isempty(info.(fieldName))
        xy = info.(fieldName);
    else
        xy = transformPoints(sourceXY, x);
    end
end

function outXY = transformPoints(sourceXY, x)
    c = cos(x(3));
    s = sin(x(3));
    outXY = sourceXY * [c, s; -s, c] + x(1:2).';
end

function value = getField(s, name, defaultValue)
    value = defaultValue;
    if isfield(s, name) && ~isempty(s.(name))
        value = s.(name);
    end
end

function angle = wrapAngle(angle)
    angle = atan2(sin(angle), cos(angle));
end
