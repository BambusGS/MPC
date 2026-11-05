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

% 5.4 Discrete-time models
disc = discretizeModels(lin, aug, p, sim.Ts);

% 5.5 Markov parameters (identical for all three models)
nMarkov = 60; % Number of Markov parameters, 60*Ts = 10 min
H = computeMarkovParameters(disc.det.A, disc.det.B, disc.det.C, nMarkov);

plotMarkovParameters("Problem 5: Markov Parameters (T_s = 10 s)", sim.Ts, H, "figures/problem5_markov_parameters.png");
