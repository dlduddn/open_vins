function frame = buildQueryBezierFrame(cfg, measFrame)
%BUILDQUERYBEZIERFRAME Convert one measurement frame to usable Bezier curves.
%
% The input pointsBody rows are Bezier control points
% [x_forward, y_left, width, curve_id]. Raw control points are never matched
% directly; every accepted item is sampled from a cubic Bezier centerline.

    frame = struct( ...
        't_cam', [], ...
        'yamlName', '', ...
        'curves', emptyCurve(), ...
        'usableIdx', [], ...
        'hasUsableCurves', false, ...
        'rejectSummary', struct(), ...
        'isYawing', false, ...
        'isPitching', false);

    if nargin < 2 || isempty(measFrame)
        return;
    end

    if isfield(measFrame, 't_cam'), frame.t_cam = measFrame.t_cam; end
    if isfield(measFrame, 'yamlName'), frame.yamlName = measFrame.yamlName; end
    if isfield(measFrame, 'isYawing'), frame.isYawing = logical(measFrame.isYawing); end
    if isfield(measFrame, 'isPitching'), frame.isPitching = logical(measFrame.isPitching); end

    if ~isfield(measFrame, 'pointsBody') || isempty(measFrame.pointsBody)
        frame.rejectSummary.empty = 1;
        return;
    end

    pointsBody = measFrame.pointsBody;
    pointsBody = pointsBody(all(isfinite(pointsBody(:, 1:min(3, size(pointsBody, 2)))), 2), :);
    if size(pointsBody, 1) < 4
        frame.rejectSummary.tooFewControlPoints = size(pointsBody, 1);
        return;
    end

    if size(pointsBody, 2) < 3
        pointsBody(:, 3) = NaN;
    end
    if size(pointsBody, 2) < 4
        pointsBody(:, 4) = 1;
    end

    curveIds = unique(pointsBody(:, 4), 'stable');
    curves = emptyCurve();
    for i = 1:numel(curveIds)
        id = curveIds(i);
        rows = pointsBody(pointsBody(:, 4) == id, :);
        curves = appendCurves(curves, controlRowsToCurves(cfg, rows, id));
    end

    if isempty(curves)
        frame.rejectSummary.noCubicBezier = 1;
        return;
    end

    curves = annotateStartConnections(cfg, curves);
    curves = applyQueryOnlyRules(cfg, curves, frame.isYawing, frame.isPitching);

    usable = find([curves.isUsable]);
    frame.curves = curves;
    frame.usableIdx = usable;
    frame.hasUsableCurves = ~isempty(usable);
    frame.rejectSummary = summarizeRejects(curves);
end

function curves = controlRowsToCurves(cfg, rows, id)
    curves = emptyCurve();
    control = rows(:, 1:2);
    widths = rows(:, 3);
    n = size(control, 1);

    if n == 4
        curves = makeBezierCurve(cfg, id, 1, control, widths);
    elseif n > 4 && mod(n - 1, 3) == 0
        segIdx = 1;
        for k = 1:3:(n - 3)
            curves = appendCurves(curves, ...
                makeBezierCurve(cfg, id, segIdx, control(k:k+3, :), widths(k:k+3)));
            segIdx = segIdx + 1;
        end
    end
end

function curve = makeBezierCurve(cfg, id, segmentIdx, controlXY, widths)
    t = linspace(0, 1, chooseSampleCount(cfg, controlXY)).';
    omt = 1 - t;
    sampleXY = ...
        (omt.^3) .* controlXY(1, :) + ...
        (3 * omt.^2 .* t) .* controlXY(2, :) + ...
        (3 * omt .* t.^2) .* controlXY(3, :) + ...
        (t.^3) .* controlXY(4, :);

    width = median(widths(isfinite(widths) & widths > 0));
    if isempty(width) || ~isfinite(width)
        width = NaN;
    end

    segLen = sqrt(sum(diff(sampleXY, 1, 1).^2, 2));
    curve = struct( ...
        'id', id, ...
        'segmentIdx', segmentIdx, ...
        'type', 'cubicBezier', ...
        'controlXY', controlXY, ...
        'sampleXY', sampleXY, ...
        'width', width, ...
        'length', sum(segLen), ...
        'chordLength', norm(controlXY(end, :) - controlXY(1, :)), ...
        'minBodyDistance', min(sqrt(sum(sampleXY.^2, 2))), ...
        'startXY', controlXY(1, :), ...
        'endXY', controlXY(end, :), ...
        'startConnected', false, ...
        'isUsable', true, ...
        'rejectReason', "");
end

function curves = annotateStartConnections(cfg, curves)
    if isempty(curves)
        return;
    end

    tol = getScalar(cfg, 'assocConnectionTol', 1.5);
    starts = vertcat(curves.startXY);
    ends = vertcat(curves.endXY);

    for i = 1:numel(curves)
        other = true(numel(curves), 1);
        other(i) = false;
        otherPts = [starts(other, :); ends(other, :)];
        if isempty(otherPts)
            curves(i).startConnected = false;
        else
            d = sqrt(sum((otherPts - curves(i).startXY).^2, 2));
            curves(i).startConnected = any(d <= tol);
        end
    end
end

function curves = applyQueryOnlyRules(cfg, curves, isYawing, isPitching)
    minChord = getScalar(cfg, 'assocMinCurveChord', 4.0);
    minBodyDist = getScalar(cfg, 'assocMinBodyDistance', 2.0);
    requireConnection = getLogical(cfg, 'assocRequireStartConnection', true);
    allowSingleCurve = getLogical(cfg, 'assocAllowSingleCurveWithoutConnection', true);
    if allowSingleCurve && numel(curves) == 1
        requireConnection = false;
    end

    for i = 1:numel(curves)
        reasons = strings(0, 1);

        if isYawing
            reasons(end + 1) = "yawing"; %#ok<AGROW>
        end
        if isPitching
            reasons(end + 1) = "pitching"; %#ok<AGROW>
        end
        if curves(i).chordLength < minChord
            reasons(end + 1) = "shortChord"; %#ok<AGROW>
        end
        if curves(i).minBodyDistance < minBodyDist
            reasons(end + 1) = "nearBody"; %#ok<AGROW>
        end
        if requireConnection && ~curves(i).startConnected
            reasons(end + 1) = "isolatedStart"; %#ok<AGROW>
        end

        curves(i).isUsable = isempty(reasons);
        if isempty(reasons)
            curves(i).rejectReason = "";
        else
            curves(i).rejectReason = strjoin(reasons, ",");
        end
    end
end

function summary = summarizeRejects(curves)
    summary = struct('usable', 0);
    if isempty(curves)
        return;
    end

    summary.usable = sum([curves.isUsable]);
    reasons = string({curves.rejectReason});
    reasons = reasons(reasons ~= "");
    for i = 1:numel(reasons)
        parts = split(reasons(i), ",");
        for j = 1:numel(parts)
            name = matlab.lang.makeValidName(char(parts(j)));
            if ~isfield(summary, name)
                summary.(name) = 0;
            end
            summary.(name) = summary.(name) + 1;
        end
    end
end

function nSamples = chooseSampleCount(cfg, controlXY)
    spacing = getScalar(cfg, 'assocCurveSampleSpacing', getScalar(cfg, 'sampleSpacing', 0.5));
    controlLen = sum(sqrt(sum(diff(controlXY, 1, 1).^2, 2)));
    nSamples = ceil(controlLen / max(spacing, eps)) + 1;
    nSamples = max(getScalar(cfg, 'assocMinCurveSamples', 12), nSamples);
    nSamples = min(getScalar(cfg, 'assocMaxCurveSamples', 80), nSamples);
end

function curves = appendCurves(curves, newCurves)
    if isempty(newCurves)
        return;
    end
    if isempty(curves)
        curves = newCurves;
    else
        curves = [curves, newCurves];
    end
end

function curve = emptyCurve()
    curve = struct( ...
        'id', {}, ...
        'segmentIdx', {}, ...
        'type', {}, ...
        'controlXY', {}, ...
        'sampleXY', {}, ...
        'width', {}, ...
        'length', {}, ...
        'chordLength', {}, ...
        'minBodyDistance', {}, ...
        'startXY', {}, ...
        'endXY', {}, ...
        'startConnected', {}, ...
        'isUsable', {}, ...
        'rejectReason', {});
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
