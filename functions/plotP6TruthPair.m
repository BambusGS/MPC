function plotP6TruthPair(t, truth, staticEstimate, dynamicEstimate, stateLabels, titleText, filePath)
    figure("Visible", "off");
    layout = tiledlayout(ceil(size(truth, 1) / 2), 2, "TileSpacing", "compact");
    for i = 1:size(truth, 1)
        nexttile;
        plot(t, truth(i, :), "k-", "LineWidth", 1.2);
        hold on;
        plot(t, staticEstimate(i, :), "--", "LineWidth", 1.0);
        plot(t, dynamicEstimate(i, :), "-", "LineWidth", 1.0);
        yline(0, ":", "Color", [0.4, 0.4, 0.4]);
        grid on;
        xlabel("Time [s]");
        ylabel(stateLabels(i));
        legend("Physical truth", "Static", "Dynamic", "Location", "best");
    end
    title(layout, titleText);
    exportgraphics(gcf, filePath, "Resolution", 400);
    close(gcf);
end
