function plotStateEstimates(t, xTrue, results, labels, stateNames, plotTitle, savePath)
% PLOTSTATEESTIMATES  True states vs Kalman filter estimates with a +-2 sigma band.
    n = size(results{1}.xf, 1);
    if nargin < 5 || isempty(stateNames)
        stateNames = arrayfun(@(i) sprintf('x_%d', i), 1:n, 'UniformOutput', false);
    end
    if nargin < 6, plotTitle = ''; end
    if nargin < 7, savePath = ''; end
    colors = [0.0000 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.4660 0.6740 0.1880; 0.4940 0.1840 0.5560];

    fig = figure('Position', [100 100 1000 220 * ceil(n / 2) + 60]);
    for i = 1:n
        subplot(ceil(n / 2), 2, i);
        hold on; box on; grid on;
        xi = results{1}.xf(i, :);
        sdi = results{1}.sdf(i, :);
        hBand = fill([t, fliplr(t)], [xi + 2 * sdi, fliplr(xi - 2 * sdi)], colors(1, :), ...
            'FaceAlpha', 0.15, 'EdgeColor', 'none');
        hEst = graphicsArray(numel(results));
        for j = 1:numel(results)
            hEst(j) = plot(t, results{j}.xf(i, :), '-', 'Color', colors(j, :), 'LineWidth', 1.0);
        end
        if ~isempty(xTrue)
            hTrue = plot(t, xTrue(i, :), 'k-', 'LineWidth', 1.2);
        end
        ylabel(stateNames{i});
        if i > n - 2, xlabel('Time'); end
        if i == 1
            if ~isempty(xTrue)
                handles = [hTrue, hEst, hBand];
                legendText = [{'True'}, labels(:)', {'\pm2\sigma'}];
            else
                handles = [hEst, hBand];
                legendText = [labels(:)', {'\pm2\sigma'}];
            end
            legend(handles, legendText, 'Location', 'best');
        end
    end
    saveFigure(fig, savePath, plotTitle);
end
