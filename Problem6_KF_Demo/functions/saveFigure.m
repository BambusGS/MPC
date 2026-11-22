function saveFigure(fig, savePath, plotTitle)
% SAVEFIGURE  Add a figure title and save the figure.
    if nargin >= 3 && ~isempty(plotTitle)
        if exist('sgtitle', 'file') || exist('sgtitle', 'builtin')
            sgtitle(fig, char(plotTitle));
        else
            annotation(fig, 'textbox', [0 0.96 1 0.04], 'String', char(plotTitle), ...
                'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
        end
    end
    if nargin < 2 || isempty(savePath)
        return
    end
    savePath = char(savePath);
    folder = fileparts(savePath);
    if ~isempty(folder) && ~exist(folder, 'dir')
        mkdir(folder);
    end
    if exist('exportgraphics', 'file') || exist('exportgraphics', 'builtin')
        exportgraphics(fig, savePath, 'Resolution', 300);
    else
        print(fig, savePath, '-dpng', '-r150');
    end
end
