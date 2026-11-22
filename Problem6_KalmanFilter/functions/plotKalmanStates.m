function plotKalmanStates(plotTitle, t, Xtrue, estimates, model, savePath)
    % True states vs Kalman filter estimates, one panel per state
    %
    % Inputs:
    %   plotTitle: Figure title
    %   t: Sample times [s], 1 x N
    %   Xtrue: n x N true states in physical units (truthInModelStates())
    %   estimates: Structure array with the fields
    %       x: n x N estimated states (model units, converted with model.scale)
    %       sd: n x N standard deviations ([] to skip the +-2 sigma band)
    %       label: Legend text (char)
    %   model: Kalman filter model with stateNames, stateUnits, scale
    %   savePath: Optional file name

    if ~isstruct(model) || ~isfield(model, 'scale')
        error('plotKalmanStates:model', 'The 5th input must be a Kalman filter model structure.');
    end
    n = size(Xtrue, 1);
    tMin = t / 60;
    colors = [0.0000 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.4660 0.6740 0.1880; 0.4940 0.1840 0.5560];

    fig = figure('Visible', 'off', 'Position', [100 100 1000 220 * ceil(n / 2) + 60]);
    for i = 1:n
        subplot(ceil(n / 2), 2, i);
        hold on; box on; grid on;

        % +-2 sigma band of the first estimate that has one
        hBand = [];
        for j = 1:numel(estimates)
            if ~isempty(estimates(j).sd)
                xi = model.scale(i) * estimates(j).x(i, :);
                sdi = model.scale(i) * estimates(j).sd(i, :);
                hBand = fill([tMin, fliplr(tMin)], [xi + 2 * sdi, fliplr(xi - 2 * sdi)], ...
                    colors(j, :), 'FaceAlpha', 0.15, 'EdgeColor', 'none');
                break
            end
        end

        hTrue = plot(tMin, Xtrue(i, :), 'k-', 'LineWidth', 1.2);
        hEst = graphicsArray(numel(estimates));
        for j = 1:numel(estimates)
            hEst(j) = plot(tMin, model.scale(i) * estimates(j).x(i, :), '-', ...
                'Color', colors(j, :), 'LineWidth', 1.0);
        end

        ylabel(sprintf('%s [%s]', model.stateNames{i}, model.stateUnits{i}));
        if i > n - 2
            xlabel('Time [min]');
        end
        if i == 1
            labels = [{'True'}, {estimates.label}];
            if isempty(hBand)
                handles = [hTrue, hEst];
            else
                handles = [hTrue, hEst, hBand];      % concatenate in one step (Line + Patch)
                labels{end + 1} = '\pm2\sigma';
            end
            legend(handles, labels, 'Location', 'northwest');
        end
    end

    saveFigure(fig, savePath, plotTitle);
    close(fig);
end
