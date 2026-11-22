function plotKalmanGains(plotTitle, t, resultDynamic, Kstatic, model, savePath)
    % Dynamic gain Kfx,k converging to the static gain Kfx
    %
    % Left: gain elements of the two measured states in physical units
    % (solid: dynamic, dashed: static). Right: ||Kfx,k - Kfx|| / ||Kfx|| (log).
    %
    % Inputs:
    %   resultDynamic: ONE kalmanFilter() result of the dynamic filter
    %   Kstatic: Static gain Kfx from kalmanStatic()
    %   model: Kalman filter model structure (uses model.scale)

    if ~isstruct(resultDynamic) || numel(resultDynamic) ~= 1 || ~isfield(resultDynamic, 'K')
        error('plotKalmanGains:result', 'The 3rd input must be one kalmanFilter() result.');
    end
    if ~isstruct(model) || ~isfield(model, 'scale')
        error('plotKalmanGains:model', 'The 5th input must be a Kalman filter model structure.');
    end

    tMin = t / 60;
    N = numel(t);
    colors = [0.0000 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.4660 0.6740 0.1880; 0.4940 0.1840 0.5560];

    fig = figure('Visible', 'off', 'Position', [100 100 1000 340]);
    subplot(1, 2, 1);
    hold on; box on; grid on;
    labels = cell(1, 4);
    c = 0;
    for i = 1:2
        for j = 1:2
            c = c + 1;
            plot(tMin, model.scale(i) * squeeze(resultDynamic.K(i, j, :)), '-', ...
                'Color', colors(c, :), 'LineWidth', 1.1);
            labels{c} = sprintf('K_{fx}(%d,%d)', i, j);
        end
    end
    c = 0;
    for i = 1:2
        for j = 1:2
            c = c + 1;
            plot(tMin, model.scale(i) * Kstatic(i, j) * ones(1, N), '--', ...
                'Color', colors(c, :), 'LineWidth', 1.0);
        end
    end
    legend(labels, 'Location', 'northeast');
    xlabel('Time [min]');
    ylabel('Gain element [cm/cm]');
    title('Dynamic (solid) and static (dashed) gain');

    subplot(1, 2, 2);
    relDiff = zeros(1, N);
    for k = 1:N
        relDiff(k) = norm(resultDynamic.K(:, :, k) - Kstatic, 'fro') / norm(Kstatic, 'fro');
    end
    semilogy(tMin, max(relDiff, 1e-16), 'LineWidth', 1.2);
    grid on;
    xlabel('Time [min]');
    ylabel('||K_{fx,k} - K_{fx}|| / ||K_{fx}||');
    title('Convergence of the dynamic gain');

    saveFigure(fig, savePath, plotTitle);
    close(fig);
end
