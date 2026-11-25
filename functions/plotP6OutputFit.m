function plotP6OutputFit(t, trueOutput, measuredOutput, C, staticResult, dynamicResult, titleText, filePath)
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
        grid on;
        xlabel("Time [min]");
        ylabel(sprintf("Y_%d [cm]", i));
        if i == 1
            legend("Truth", "Measured", "Static", "Dynamic", "Location", "best");
        end
    end
    title(layout, titleText);
    exportgraphics(gcf, filePath, "Resolution", 400);
    close(gcf);
end
