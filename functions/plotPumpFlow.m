function plotPumpFlow(plotTitle, time, pumpFlows, savePath, flowNames)
    % Plot the two pump flow rates.
    %
    % pumpFlows must have one row per flow signal.
    % flowNames and savePath are optional.

    figure;
    clf;

    sgtitle(plotTitle);

    % MATLAB default colors
    colors = [
        0.0000 0.4470 0.7410;
        0.8500 0.3250 0.0980
    ];

    for pumpIndex = 1:2
        subplot(2, 1, pumpIndex);
        plot(time / 60, pumpFlows(pumpIndex, :).', "Color", colors(pumpIndex, :), "LineWidth", 1.5);
        xlabel("Time [min]");
        if nargin >= 5 && ~isempty(flowNames)
            flowName = flowNames{pumpIndex};
        else
            flowName = sprintf("Flow %d", pumpIndex);
        end
        ylabel(sprintf("%s [cm^3/s]", flowName));
        grid on;
    end

    if nargin >= 4 && ~isempty(savePath)
        exportgraphics(gcf, savePath, "Resolution", 400);
    end
end