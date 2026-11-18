%% Setup
% Problem 5 uses the identified models of Problem 4 (step_sizes, estimated_params, G_final).
% Run problem4_main.m first; this script does not clear the workspace.

addpath(fullfile(fileparts(mfilename("fullpath")), "functions"));
% Standalone script for Problem 5, to be merged into main.m later.

p = createParameters();
sim = createSimSettings(p);

%% Problem 5

% 5.1 Steady state and continuous-time linearized models
op = computeSteadyState([sim.F1; sim.F2], p.FMean, p); % Operating point
lin = computeJacobians(op, p);                         % Deterministic and piecewise model
aug = computeAugmentedModel(lin, p);                   % SDE model with disturbance states

% 5.2 Continuous-time transfer functions
trf = computeTransferFunctions(lin, p);

% 5.3 Comparison with the models identified from step responses (Problem 4)
% estimated_params(tank, pump, :, step) = [K, Tmain, Textra], with Textra = 0 on the diagonal
ident.stepSizes = step_sizes;
ident.K = squeeze(estimated_params(:, :, 1, :));
ident.T = cell(1, numel(step_sizes));
for s = 1:numel(step_sizes)
    ident.T{s} = cell(2);
    for i = 1:2
        for j = 1:2
            T = squeeze(estimated_params(i, j, 2:3, s)).';
            ident.T{s}{i, j} = T(T > 0);
        end
    end
end

cmp = compareIdentifiedModel(trf, ident);
disp(cmp);

% 5.4 Discrete-time models
disc = discretizeModels(lin, aug, p, sim.Ts);

% 5.5 Markov parameters (identical for all three models)
nMarkov = 90; % 90*Ts = 15 min
H = computeMarkovParameters(disc.det.A, disc.det.B, disc.det.C, nMarkov);

% Markov parameters of the identified models from Problem 4
Hid = cell(1, numel(G_final));
for s = 1:numel(G_final)
    [A, B, C] = ssdata(c2d(ss(G_final{s}), sim.Ts, "zoh"));
    Hid{s} = computeMarkovParameters(A, B, C, nMarkov);
end

plotMarkovParameters("Markov Parameters (T_s = 10 s)", sim.Ts, [{H}, Hid], "figures/problem5_markov_parameters.png", ...
    ["Linearized", "Identified 10%", "Identified 25%", "Identified 50%"]);
