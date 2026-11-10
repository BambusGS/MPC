function cmp = compareIdentifiedModel(trf, ident)
    % Compare gains and time constants of the linearized model G(s) with
    % those identified from the step responses (Problem 4)
    %
    % Inputs:
    %   trf: Structure created by computeTransferFunctions() (fields Ku, Tu)
    %   ident: Structure with the identified models
    %       stepSizes: 1 x nSteps relative step sizes [-]
    %       K: 2 x 2 x nSteps gains [cm/(cm^3/s)]
    %       T: 1 x nSteps cell, each a 2x2 cell of time constants [s],
    %          ordered as in trf.Tu ([bottom tank, top tank] off the diagonal)
    %
    % Output:
    %   cmp: Table with one row per parameter, the linearized value and the
    %        identified value and relative deviation [%] for each step size

    nSteps = numel(ident.stepSizes);
    names = strings(0, 1);
    lin = zeros(0, 1);
    idv = zeros(0, nSteps);

    for i = 1:2         % Output (tank)
        for j = 1:2     % Input (pump)
            % Gain
            names(end + 1, 1) = sprintf("K%d%d", i, j);
            lin(end + 1, 1) = trf.Ku(i, j);
            idv(end + 1, :) = squeeze(ident.K(i, j, :)).';

            % Time constants
            for l = 1:numel(trf.Tu{i, j})
                names(end + 1, 1) = sprintf("T%d%d_%d", i, j, l);
                lin(end + 1, 1) = trf.Tu{i, j}(l);
                idv(end + 1, :) = cellfun(@(T) T{i, j}(l), ident.T);
            end
        end
    end

    dev = 100 * (idv - lin) ./ lin; % Relative deviation [%]

    % Assemble table: Lin | Id10 | Dev10 | Id25 | Dev25 | ...
    cmp = table(names, lin, 'VariableNames', ["Parameter", "Linearized"]);
    for s = 1:nSteps
        pct = round(100 * ident.stepSizes(s));
        cmp.(sprintf("Id%d", pct)) = idv(:, s);
        cmp.(sprintf("Dev%d", pct)) = dev(:, s);
    end
end
