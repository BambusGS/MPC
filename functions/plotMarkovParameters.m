function plotMarkovParameters(plotTitle, Ts, H, savePath, labels)

    % Plot the Markov parameters of a 2x2 system, one subplot per input-output pair.
    %
    % H is a 2-by-2-by-(N+1) array from computeMarkovParameters(), or a cell
    % array of such arrays to compare several models (first one as stems).

    if ~iscell(H)
        H = {H};
    end

    figure;
    sgtitle(plotTitle);
    time = (0:size(H{1}, 3) - 1) * Ts / 60;

    for i = 1:2
        for j = 1:2
            subplot(2, 2, 2 * (i - 1) + j);
            hold on;
            stem(time, squeeze(H{1}(i, j, :)), "filled", "MarkerSize", 3);
            for m = 2:numel(H)
                plot(time, squeeze(H{m}(i, j, :)), "LineWidth", 1.2);
            end
            hold off;
            title(sprintf("F%d \\rightarrow h%d", j, i));
            xlabel("Time [min]");
            ylabel("H_i [cm/(cm^3/s)]");
            grid on;
        end
    end

    if nargin >= 5
        subplot(2, 2, 2);
        legend(labels);
    end
    if nargin >= 4 && ~isempty(savePath)
        exportgraphics(gcf, savePath, "Resolution", 400);
    end
end
