%% Setup
clear;
clc;
close all;

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
% TODO: take these from the Problem 4 results struct once merged into main.m
ident.stepSizes = [0.10 0.25 0.50];
ident.K = cat(3, [0.133583 0.065897;  0.118061 0.20865],  ...
                 [0.139018 0.067275;  0.120523 0.216012], ...
                 [0.148062 0.0695735; 0.124635 0.228254]);
ident.T = {{85.6347, [85.4530 35.9792]; [115.362 56.0529], 114.452}, ...
           {90.6103, [90.1486 35.9275]; [122.356 55.3614], 120.191}, ...
           {98.9596, [98.0254 35.8703]; [133.882 54.4624], 129.790}};

cmp = compareIdentifiedModel(trf, ident);
disp(cmp);

% 5.4 Discrete-time models
disc = discretizeModels(lin, aug, p, sim.Ts);

% 5.5 Markov parameters (identical for all three models)
nMarkov = 90; % 90*Ts = 15 min
H = computeMarkovParameters(disc.det.A, disc.det.B, disc.det.C, nMarkov);

% Markov parameters of the identified models from Problem 4
% TODO: take G_final from Problem 4 once merged into main.m
Hid = cell(1, 3);
for s = 1:3
    G = tf(zeros(2));
    for i = 1:2
        for j = 1:2
            den = 1;
            for T = ident.T{s}{i, j}
                den = conv(den, [T 1]);
            end
            G(i, j) = tf(ident.K(i, j, s), den);
        end
    end
    [A, B, C] = ssdata(c2d(ss(G), sim.Ts, "zoh"));
    Hid{s} = computeMarkovParameters(A, B, C, nMarkov);
end

plotMarkovParameters("Markov Parameters (T_s = 10 s)", sim.Ts, [{H}, Hid], "figures/problem5_markov_parameters.png", ...
    ["Linearized", "Identified 10%", "Identified 25%", "Identified 50%"]);
