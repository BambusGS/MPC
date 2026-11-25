function plotP6OutputFit(t, trueOutput, measuredOutput, C, staticResult, dynamicResult, titleText, filePath)
    staticPrediction = C * staticResult.xp;
    dynamicPrediction = C * dynamicResult.xp;
    figure("Visible", "off");
    layout = tiledlayout(size(C, 1), 1, "TileSpacing", "compact");
    for i = 1:size(C, 1)
        nexttile;
        plot(t, trueOutput(i, :), "k-", "LineWidth", 1.2);
        hold on;
        plot(t, measuredOutput(i, :), ".", "Color", [0.65, 0.65, 0.65], "MarkerSize", 5);
        plot(t, staticPrediction(i, :), "--", "LineWidth", 1.0);
        plot(t, dynamicPrediction(i, :), "-", "LineWidth", 1.0);
        grid on;
        xlabel("Time [s]");
        ylabel(sprintf("Y_%d [cm]", i));
        legend("True output", "Measured", "Static prediction", "Dynamic prediction", ...
            "Location", "best");
    end
    title(layout, titleText);
    exportgraphics(gcf, filePath, "Resolution", 400);
    close(gcf);
end
