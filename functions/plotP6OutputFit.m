function plotP6OutputFit(t, trueOutput, measuredOutput, C, staticResult, dynamicResult, titleText, truthLabel, filePath, stepTimes)
    if nargin < 10
        stepTimes = [];
    end

    staticPrediction = C * staticResult.xp;
    dynamicPrediction = C * dynamicResult.xp;
    figure("Visible", "off");
    layout = tiledlayout(size(C, 1), 1, "TileSpacing", "compact");
    for i = 1:size(C, 1)
        nexttile;
        plot(t / 60, trueOutput(i, :), "-", "LineWidth", 1.2);
        hold on;
        plot(t / 60, measuredOutput(i, :), ".", "MarkerSize", 5);
        plot(t / 60, staticPrediction(i, :), "--", "LineWidth", 1.0);
        plot(t / 60, dynamicPrediction(i, :), "-", "LineWidth", 1.0);
        for stepIndex = 1:numel(stepTimes)
            xline(stepTimes(stepIndex) / 60, "--", "Color", "k", "HandleVisibility", "off");
        end
        grid on;
        xlabel("Time [min]");
        ylabel(sprintf("Y_%d [cm]", i));
        if i == 1
            legend(truthLabel, "Measured", "Static KF", "Dynamic KF", "Location", "best");
        end
    end
    title(layout, titleText);
    exportgraphics(gcf, filePath, "Resolution", 400);
    close(gcf);
end
