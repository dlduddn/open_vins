function metrics = evaluateTrajectories(gtXY, vioXY, pfXY)
%EVALUATETRAJECTORIES Compare VIO and PF positions over the same trajectory.

    n = min([size(gtXY, 2), size(vioXY, 2), size(pfXY, 2)]);
    gtXY = gtXY(:, 1:n);
    vioXY = vioXY(:, 1:n);
    pfXY = pfXY(:, 1:n);

    metrics = struct();
    metrics.n = n;
    metrics.vio = positionMetrics(vioXY, gtXY);
    metrics.pf = positionMetrics(pfXY, gtXY);
    metrics.rmseImprovement = metrics.vio.rmse - metrics.pf.rmse;
    metrics.rmseImprovementPct = 100.0 * metrics.rmseImprovement / max(metrics.vio.rmse, eps);
end

function out = positionMetrics(estXY, gtXY)
    valid = all(isfinite(estXY), 1) & all(isfinite(gtXY), 1);
    err = estXY(:, valid) - gtXY(:, valid);
    euclidean = sqrt(sum(err .^ 2, 1));

    out = struct();
    out.validCount = nnz(valid);
    if isempty(euclidean)
        out.rmse = NaN;
        out.mean = NaN;
        out.median = NaN;
        out.final = NaN;
        out.axisRmse = [NaN; NaN];
        out.error = euclidean;
        return;
    end

    out.rmse = sqrt(mean(euclidean .^ 2, 'omitnan'));
    out.mean = mean(euclidean, 'omitnan');
    out.median = median(euclidean, 'omitnan');
    out.final = euclidean(end);
    out.axisRmse = sqrt(mean(err .^ 2, 2, 'omitnan'));
    out.error = euclidean;
end
