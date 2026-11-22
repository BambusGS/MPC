function plotInnovations(plotTitle, t, results, labels, savePath)
    % Innovations of several Kalman filters with +-2*sqrt(Re) of the first one.
    % A nonzero mean after a disturbance step shows a steady-state offset.
    %
    % Inputs:
    %   results: Cell array of kalmanFilter() results
    %   labels: Cell array of legend texts (char)

    tMin = t / 60;
    colors = [0.0000 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.4660 0.6740 0.1880; 0.4940 0.1840 0.5560];

    fig = figure('Visible', 'off', 'Position', [100 100 1000 430]);
    for i = 1:2
        subplot(2, 1, i);
        hold on; box on; grid on;
        h = graphicsArray(numel(results));
        for j = 1:numel(results)
            h(j) = plot(tMin, results{j}.e(i, :), '-', 'Color', colors(j, :), 'LineWidth', 0.9);
        end
        sd = sqrt(squeeze(results{1}.Re(i, i, :)))';
        plot(tMin, 2 * sd, 'k:', tMin, -2 * sd, 'k:');
        ylabel(sprintf('e_%d [cm]', i));
        if i == 1
            legend(h, labels, 'Location', 'northwest');
        end
    end
    xlabel('Time [min]');

    saveFigure(fig, savePath, plotTitle);
    close(fig);
end
