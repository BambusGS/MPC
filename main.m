%% Setup
clear;
clc;
close all;

addpath(fullfile(fileparts(mfilename("fullpath")), "functions"));
% all functions live in the "functions" folder. Pass data in and return data out. Avoid relying on global variables.
% for some reason you need to run this script with F5, otherwise it does
%   not find the functions/ folder ¯\_(ツ)_/¯

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

%% Problem 7, simple known quadratic program (QP) solver test

% Solve: Quadratic Minimization with Linear Constraints and Bounds [MATLAB example]
%
%   min  phi = 1/2*x'*H*x + g'*x
%    x
%
%   H = [ 1 -1  1
%        -1  2 -2
%         1 -2  4]
%
%   g = [2; -3; 1]
%
% Subject to:
%
%   0 <= x <= 1
%   x1 + x2 + x3 = 1/2
%
% The equality constraint sum(x) = 1/2 is represented by:
%   1/2 <= [1 1 1]*x <= 1/2    so that it fits the required solver form

H = [1 -1  1;
    -1  2 -2;
     1 -2  4];

g = [2; -3; 1];

l = zeros(3, 1);
u = ones(3, 1);

A = ones(1, 3);
bl = 1/2;
bu = 1/2;

xinit = [];

[x, info] = qpsolver(H, g, l, u, A, bl, bu, xinit); % internal options exit

expectedX = [0; 0.5; 0]; % from MATLAB site
disp("QP solution:");
disp(x);

disp("QP solver info:");
disp(info);

assert(info.exitflag > 0, "QP solver failed.");
assert(norm(x - expectedX) < 1e-9, "QP solution is incorrect.");

disp("QP test passed.");

%% Problem 8