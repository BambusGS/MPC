function plotMarkovParameters(plotTitle, Ts, H, savePath)

    % Plot the Markov parameters of a 2x2 system, one subplot per input-output pair.
    %
    % H must be 2-by-2-by-(N+1) with H(:, :, i+1) = H_i, as returned by
    % computeMarkovParameters().
    %
    % savePath is optional.

    figure;
    clf;

    sgtitle(plotTitle);

    time = (0:size(H, 3) - 1) * Ts;

    for outputIndex = 1:2
        for inputIndex = 1:2
            subplot(2, 2, 2 * (outputIndex - 1) + inputIndex);
            stem(time / 60, squeeze(H(outputIndex, inputIndex, :)), "filled", "MarkerSize", 3, "LineWidth", 1.2);

            title(sprintf("F%d \\rightarrow h%d", inputIndex, outputIndex));
            xlabel("Time [min]");
            ylabel("H_i [cm/(cm^3/s)]");
            grid on;
        end
    end

    if nargin >= 4 && ~isempty(savePath)
        exportgraphics(gcf, savePath, "Resolution", 400);
    end
end
