function stats = innovationStats(result, idx)
    % Checks of a Kalman filter through its innovations
    %
    % For a correct model the innovations are white, zero-mean and have the
    % covariance Re (lecture slide 38, innovation form).
    %
    % Inputs:
    %   result: Output of kalmanFilter()
    %   idx: Sample indices to use (e.g. after the initial transient)
    %
    % Output:
    %   stats: mean, cov (sample covariance), Re (mean predicted covariance),
    %       NIS (mean normalized innovation squared, expected value = ny),
    %       rho1 (lag-1 autocorrelation, approx. 0 if white)

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
