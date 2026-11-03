%% Setup
clear;
clc;
close all;

addpath(fullfile(fileparts(mfilename("fullpath")), "functions"));
% all functions live in the "functions" folder. Pass data in and return data out. Avoid relying on global variables.

p = createParameters();
sim = createSimSettings(p);

%% Problem 1

%% Problem 2

rng(sim.randomSeed);

% Simulate the three nonlinear open-loop model variants.
deterministic = simulateDeterministic(p, sim);

rng(sim.randomSeed);
piecewise = simulatePiecewise(p, sim);

rng(sim.randomSeed);
sde = simulateSDE(p, sim);

% Full state is used only for diagnostic tank-height plots.
deterministicHeights = deterministic.x ./ (p.rho .* p.A.');
piecewiseHeights = piecewise.x ./ (p.rho .* p.A.');
sdeHeights = sde.x(:, 1:4) ./ (p.rho .* p.A.');


plotTankHeights("Problem 2: Constant Disturbances - Tank Heights", deterministic.t, deterministicHeights, "figures/problem2_deterministic_tanks.png");
plotPumpFlow("Problem 2: Constant Disturbances - Pump Inputs", sim.t, sim.input, "figures/problem2_deterministic_inputs.png", {"F1", "F2"});
plotPumpFlow("Problem 2: Constant Disturbances - Disturbance Flows", deterministic.t, deterministic.d.', "figures/problem2_deterministic_disturbances.png", {"F3", "F4"});

plotTankHeights("Problem 2: Piecewise-Constant Disturbances - Tank Heights", piecewise.t, piecewiseHeights, "figures/problem2_piecewise_tanks.png");
plotPumpFlow("Problem 2: Piecewise-Constant Disturbances - Pump Inputs", sim.t, sim.input, "figures/problem2_piecewise_inputs.png", {"F1", "F2"});
plotPumpFlow("Problem 2: Piecewise-Constant Disturbances - Disturbance Flows", piecewise.tPlot, piecewise.dPlot.', "figures/problem2_piecewise_disturbances.png", {"F3", "F4"});

plotTankHeights("Problem 2: SDE Disturbances - Tank Heights", sde.t, sdeHeights, "figures/problem2_sde_tanks.png");
plotPumpFlow("Problem 2: SDE Disturbances - Pump Inputs", sim.t, sim.input, "figures/problem2_sde_inputs.png", {"F1", "F2"});
plotPumpFlow("Problem 2: SDE Disturbances - Disturbance Flows", sde.t, sde.d.', "figures/problem2_sde_disturbances.png", {"F3", "F4"});