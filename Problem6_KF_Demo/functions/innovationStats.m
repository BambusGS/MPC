function stats = innovationStats(result, idx)
% INNOVATIONSTATS  Innovation mean, covariance, NIS and lag-1 autocorrelation.
    e = result.e(:, idx);
    stats.mean = mean(e, 2);
    stats.cov = cov(e');
    stats.Re = mean(result.Re(:, :, idx), 3);

    nis = zeros(1, numel(idx));
    for k = 1:numel(idx)
        nis(k) = e(:, k)' / result.Re(:, :, idx(k)) * e(:, k);
    end
    stats.NIS = mean(nis);

    stats.rho1 = zeros(size(e, 1), 1);
    for i = 1:size(e, 1)
        ei = e(i, :) - mean(e(i, :));
        stats.rho1(i) = sum(ei(1:end-1) .* ei(2:end)) / sum(ei .^ 2);
    end
end
