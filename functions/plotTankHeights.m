function plotTankHeights(plotTitle, calculatedTime, calculatedHeights, savePath)

    % Plot the heights of all four tanks.
    %
    % calculatedHeights must have one row per time sample
    % and four columns, one column per tank.
    %
    % savePath is optional. Heights must be N-by-4.

    figure;
    clf;

    sgtitle(plotTitle);

    % Plot tanks in the order 3, 4, 1, 2
    tankOrder = [3 4 1 2];

    % MATLAB default colors
    colors = [
        0.0000 0.4470 0.7410;
        0.8500 0.3250 0.0980;
        0.9290 0.6940 0.1250;
        0.4940 0.1840 0.5560
    ];

    for plotIndex = 1:4
        tankIndex = tankOrder(plotIndex);
        subplot(2, 2, plotIndex);
        plot(calculatedTime / 60, calculatedHeights(:, tankIndex), "-", "Color", colors(tankIndex, :), "LineWidth", 1.8);

        xlabel("Time [min]");
        ylabel(sprintf("Tank %d Height [cm]", tankIndex));
        grid on;
    end

    if nargin >= 4 && ~isempty(savePath)
        exportgraphics(gcf, savePath, "Resolution", 400);
    end
end