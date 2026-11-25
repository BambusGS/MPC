function plotP6TruthPair(t, truth, staticEstimate, dynamicEstimate, stateLabels, titleText, truthLabel, filePath, stepTimes)
    if nargin < 9
        stepTimes = [];
    end

    figure("Visible", "off");
    layout = tiledlayout(ceil(size(truth, 1) / 2), 2, "TileSpacing", "compact");
    for i = 1:size(truth, 1)
        nexttile;
        plot(t / 60, truth(i, :), "-", "LineWidth", 1.2);
        hold on;
        plot(t / 60, staticEstimate(i, :), "--", "LineWidth", 1.0);
        plot(t / 60, dynamicEstimate(i, :), "-", "LineWidth", 1.0);
        for stepIndex = 1:numel(stepTimes)
            xline(stepTimes(stepIndex) / 60, "--", "Color", "k", "HandleVisibility", "off");
        end
        yline(0, ":", "Color", [0.4, 0.4, 0.4]);
        grid on;
        xlabel("Time [min]");
        ylabel(stateLabels(i));
        if i == 1
            legend(truthLabel, "Static KF", "Dynamic KF", "Location", "best");
        end
    end
    title(layout, titleText);
    exportgraphics(gcf, filePath, "Resolution", 400);
    close(gcf);
end
