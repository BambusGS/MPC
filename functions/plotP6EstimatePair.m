function plotP6EstimatePair(t, staticResult, dynamicResult, stateIndices, stateLabels, titleText, filePath)
    figure("Visible", "off");
    layout = tiledlayout(ceil(numel(stateIndices) / 2), 2, "TileSpacing", "compact");
    for i = 1:numel(stateIndices)
        nexttile;
        stateIndex = stateIndices(i);
        plot(t, staticResult.xf(stateIndex, :), "--", "LineWidth", 1.1);
        hold on;
        plot(t, dynamicResult.xf(stateIndex, :), "-", "LineWidth", 1.1);
        yline(0, ":", "Color", [0.4, 0.4, 0.4]);
        grid on;
        xlabel("Time [s]");
        ylabel(stateLabels(i));
        legend("Static", "Dynamic", "Location", "best");
    end
    title(layout, titleText);
    exportgraphics(gcf, filePath, "Resolution", 400);
    close(gcf);
end
