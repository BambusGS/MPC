%% Problem 6 setup
% Load shared physical parameters and Problem 6 helper functions.
projectRoot = fileparts(mfilename("fullpath"));
addpath(fullfile(projectRoot, "functions"));
addpath(fullfile(projectRoot, "Problem6_KF_Demo", "functions"));
p = createParameters();
sim = createSimSettings(p);

%% 6.1 Configuration, operating point, and repeatable simulations
% Define the nominal disturbance, sensor noise, and steady-state initialization.
p6 = p;
p6Results = struct();
p6.FMean = [50; 100];                     % Nominal [F3; F4] [cm^3/s]
p6.disturbanceRate = [0.05; 0.05];         % Mean-reversion rates [1/s]
p6.disturbanceSigma = [5; 5];              % Diffusion coefficients
p6.measurementCovariance = eye(2);         % Sensor covariance [cm^2]
p6.measurementVariance = 1;

% Configure a Problem 6-specific time grid without changing the shared simulation settings.
sim6 = sim;
sim6.Ts = 10;
sim6.tf = 30 * 60;
sim6.t = (sim6.t0:sim6.Ts:sim6.tf).';
sim6.N = numel(sim6.t);
operpt6 = computeSteadyState([300; 300], p6.FMean, p6);
sim6.initialState = operpt6.x;
sim6.initialDisturbance = p6.FMean;
inputDeviation6 = zeros(2, sim6.N);        % Pumps stay at their nominal flows
measurementCovariance6 = eye(2);
P0Scale6 = 5;                              % Same documented demo initialization, P0 = 5 I

% Verify that the SDE integration grid aligns with the discrete filter samples.
assert(abs(sim6.Ts / sim6.dt - round(sim6.Ts / sim6.dt)) < 1e-12, "Problem 6 sampling interval must be an integer number of SDE steps.");
sampleStride6 = round(sim6.Ts / sim6.dt);
sampleIndex6 = 1:sampleStride6:round((sim6.tf - sim6.t0) / sim6.dt) + 1;
assert(numel(sampleIndex6) == sim6.N, "SDE and KF time grids do not align.");
sampleTimes6 = sim6.t;

fprintf("Finished P6.1 Setup\n");
%% 6.2 P4 10% identified model and four-state realization
% Set the identified transfer gains and shared-pole time constants.
% The gains are the 10% identified values reported in Chapter 5. The shared time constants are the 10% Problem 4 values reported in Chapter 4.
K11 = 0.133583;
K12 = 0.0658997;
K21 = 0.118038;
K22 = 0.208652;
T1 = 85.6311;
T2 = 114.467;
T3 = 35.8382;
T4 = 56.7755;

% Continuous realization of the Chapter 4 shared-pole transfer matrix.
A4c = diag([-1 / T1, -1 / T2, -1 / T3, -1 / T4]);
A4c(1, 3) = K12 / T1;
A4c(2, 4) = K21 / T2;
B4c = [K11 / T1, 0;
       0, K22 / T2;
       0, 1 / T3;
       1 / T4, 0];
C4c = [1, 0, 0, 0;
       0, 1, 0, 0];
% Discretize the continuous P4 realization at the Problem 6 sampling interval.
[A4z, B4z] = c2dzoh(A4c, B4c, sim6.Ts);


% Generate Markov parameters and assemble the Hankel matrix and its shift.
hankelSizeP4 = 44;
nMarkovP4 = 2 * hankelSizeP4;
H4transfer = computeMarkovParameters(A4z, B4z, C4c, nMarkovP4);

HankelP4 = zeros(2 * hankelSizeP4);
HankelShiftP4 = zeros(size(HankelP4));

for blockRow = 1:hankelSizeP4
    for blockColumn = 1:hankelSizeP4
        rows = (2 * blockRow - 1):(2 * blockRow);
        columns = (2 * blockColumn - 1):(2 * blockColumn);

        % The helper stores H0 at index 1, so Hk is at index k+1.
        markovIndex = blockRow + blockColumn - 1;
        HankelP4(rows, columns) = H4transfer(:, :, markovIndex + 1);
        HankelShiftP4(rows, columns) = H4transfer(:, :, markovIndex + 2);
    end
end

[U4, S4, V4] = svd(HankelP4, "econ");

% Retain the four dominant Hankel directions.
r = 4;
Ur = U4(:, 1:r);
Vr = V4(:, 1:r);
Sr = diag(sqrt(diag(S4(1:r, 1:r))));

% Recover the minimal P4 realization from the shifted Hankel matrix.
A4id = diag(1 ./ diag(Sr)) * Ur' * HankelShiftP4 * Vr * diag(1 ./ diag(Sr));
B4id = Sr * Vr(1:2, :)';
C4id = Ur(1:2, :) * Sr;

% Check that the identified and discretized P4 models have the same poles.
assert(norm(sort(eig(A4id)) - sort(eig(A4z)), inf) < 1e-8, "P4 Hankel realization has different poles from the discretized model.");


% Define the standard discrete-time Kalman models for P4.
modelP4.A = A4id;
modelP4.B = B4id;
modelP4.G = eye(4);
modelP4.C = C4id;
modelP4.Q = 5^2 * eye(4);
modelP4.R = measurementCovariance6;
modelP4.S = zeros(4, 2);
modelP4.P0 = P0Scale6 * eye(4);
modelP4.name = "P4 identified four-state";
assertKalmanModelP6(modelP4);

% Mark the successful identification and realization checks for P4.
fprintf("Finished P6.2 P4 identified model and four-state realization\n");

%% 6.3 P5 numerical six-state SDE model
% These are Chapter 5 equation (Qbar_numbers) values, state order
% xa = [X1, X2, X3, X4, D3, D4]' with D = F - [50; 100].
A6bar = [0.886, 0,     0.228, 0,     1.032, 0;
         0,     0.914, 0,     0.155, 0,     0.687;
         0,     0,     0.758, 0,     6.793, 0;
         0,     0,     0,     0.838, 0,     7.161;
         0,     0,     0,     0,     0.607, 0;
         0,     0,     0,     0,     0,     0.607];
B6bar = [5.650, 0.365;
         0.323, 6.693;
         0,     2.620;
         3.667, 0;
         0,     0;
         0,     0];
Q6bar = [58.9, 0,    505.1, 0,    64.4, 0;
         0,    25.7, 0,     346.8, 0,    42.3;
         505.1, 0,   4717.1, 0,   707.7, 0;
         0,    346.8, 0,    5086.5, 0,   730.8;
         64.4, 0,   707.7, 0,    158.0, 0;
         0,    42.3, 0,    730.8, 0,    158.0];

% Define the standard discrete-time Kalman models for P5.
modelP5.A = A6bar;
modelP5.B = B6bar;
modelP5.G = eye(6);
modelP5.C = [lin6.C, zeros(2, 2)];
modelP5.Q = Q6bar;
modelP5.R = measurementCovariance6;
modelP5.S = zeros(6, 2);
modelP5.P0 = P0Scale6 * eye(6);
modelP5.name = "P5 physical six-state SDE";
assertKalmanModelP6(modelP5);

% Mark completion of the P4/P5 model definitions and validations.
fprintf("Finished P6.3 P5 numerical six-state SDE model\n");

%% 6.4 Problem 4 linear model, no step
% Simulate the identified model and compare static and dynamic filters on shared observations.
% Simulate with the same Q and R used to design the filters.
rng(sim.randomSeed + 42);
modelP4Truth = modelP4;
modelP4Truth.P0 = zeros(4);
truthP4LinearNoStep = simulateLinearSystem(modelP4Truth, inputDeviation6, sim6.N);
resultP4LinearNoStep.dynamic = kalmanFilter(modelP4, truthP4LinearNoStep.y, inputDeviation6, "dynamic");
resultP4LinearNoStep.static = kalmanFilter(modelP4, truthP4LinearNoStep.y, inputDeviation6, "static");
p6Results.p4.linear.noStep.truth = truthP4LinearNoStep;
p6Results.p4.linear.noStep.filters = resultP4LinearNoStep;
p6Results.p4.linear.noStep.outputRMSE = [sqrt(mean((modelP4.C * resultP4LinearNoStep.static.xf - modelP4.C * truthP4LinearNoStep.x).^2, "all")), sqrt(mean((modelP4.C * resultP4LinearNoStep.dynamic.xf - modelP4.C * truthP4LinearNoStep.x).^2, "all"))];
plotP6OutputFit(sampleTimes6, modelP4.C * truthP4LinearNoStep.x, truthP4LinearNoStep.y, modelP4.C, resultP4LinearNoStep.static, resultP4LinearNoStep.dynamic, "P4 linear no-step", fullfile(projectRoot, "figures", "problem6_4_p4_linear_no_step_outputs.png"));

fprintf("Finished P6.4 P4 linear model, no step\n");
%% 6.5 Problem 5 linear model, no step
% Simulate the physical model and compare both filters against its linear-model truth.
rng(sim.randomSeed + 43);
modelP5Truth = modelP5;
modelP5Truth.P0 = zeros(6);
truthP5LinearNoStep = simulateLinearSystem(modelP5Truth, inputDeviation6, sim6.N);
resultP5LinearNoStep.dynamic = kalmanFilter(modelP5, truthP5LinearNoStep.y, inputDeviation6, "dynamic");
resultP5LinearNoStep.static = kalmanFilter(modelP5, truthP5LinearNoStep.y, inputDeviation6, "static");
p6Results.p5.linear.noStep.truth = truthP5LinearNoStep;
p6Results.p5.linear.noStep.filters = resultP5LinearNoStep;
p6Results.p5.linear.noStep.outputRMSE = [sqrt(mean((modelP5.C * resultP5LinearNoStep.static.xf - modelP5.C * truthP5LinearNoStep.x).^2, "all")), sqrt(mean((modelP5.C * resultP5LinearNoStep.dynamic.xf - modelP5.C * truthP5LinearNoStep.x).^2, "all"))];
plotP6OutputFit(sampleTimes6, modelP5.C * truthP5LinearNoStep.x, truthP5LinearNoStep.y, modelP5.C, resultP5LinearNoStep.static, resultP5LinearNoStep.dynamic, "P5 linear no-step", fullfile(projectRoot, "figures", "problem6_5_p5_linear_no_step_outputs.png"));

fprintf("Finished P6.5 P5 linear model, no step\n");
%% 6.6 Nonlinear model, no step
% Generate one nonlinear truth realization and sample it on the filter grid.
rng(sim.randomSeed);
truthNoStepFine = simulateSDE_rerun(p6, sim6);
truthNoStep6 = truthNoStepFine.x(sampleIndex6, :);
truthNoStepD6 = truthNoStepFine.d(sampleIndex6, :);
measurementNoStep6 = truthNoStepFine.y(sampleIndex6, :).'- operpt6.y;
trueTankNoStep6 = tankHeightsFromMasses(truthNoStep6(:, 1:4).', p6) - operpt6.h;
trueDisturbanceNoStep6 = truthNoStepD6.' - p6.FMean;
trueOutputNoStep6 = truthNoStep6(:, 1:2).' ./ (p6.rho .* p6.A(1:2)) - operpt6.y;
resultP4NoStep.dynamic = kalmanFilter(modelP4, measurementNoStep6, inputDeviation6, "dynamic");
resultP4NoStep.static = kalmanFilter(modelP4, measurementNoStep6, inputDeviation6, "static");
resultP5NoStep.dynamic = kalmanFilter(modelP5, measurementNoStep6, inputDeviation6, "dynamic");
resultP5NoStep.static = kalmanFilter(modelP5, measurementNoStep6, inputDeviation6, "static");
p6Results.nonlinear.noStep.truth = truthNoStepFine;
p6Results.nonlinear.noStep.filters.p4 = resultP4NoStep;
p6Results.nonlinear.noStep.filters.p5 = resultP5NoStep;

% Confirm that both filters use the same observations and return finite states.
assert(isequal(measurementNoStep6, truthNoStepFine.y(sampleIndex6, :).'-operpt6.y) && isequal(inputDeviation6, zeros(size(inputDeviation6))), "The filters must use the same sampled measurements and known input deviations.");
assert(all(isfinite(resultP4NoStep.dynamic.xf(:))) && all(isfinite(resultP5NoStep.dynamic.xf(:))), "No-step filter returned nonfinite state estimates.");

% Measure convergence of the time-varying filter gains to their steady-state values.
[~, staticGainP4] = kalmanStatic(modelP4);
[~, staticGainP5] = kalmanStatic(modelP5);
gainGapP4 = zeros(1, sim6.N);
gainGapP5 = zeros(1, sim6.N);
for k = 1:sim6.N
    gainGapP4(k) = norm(resultP4NoStep.dynamic.K(:, :, k) - staticGainP4, "fro");
    gainGapP5(k) = norm(resultP5NoStep.dynamic.K(:, :, k) - staticGainP5, "fro");
end
gainRelativeErrorP4 = gainGapP4(end) / max(norm(staticGainP4, "fro"), eps);
gainRelativeErrorP5 = gainGapP5(end) / max(norm(staticGainP5, "fro"), eps);
assert(gainRelativeErrorP4 < 1e-5 && gainRelativeErrorP5 < 1e-5, "Dynamic no-step filter gains did not converge to their static gains.");

% Report no-step state and output errors in their corresponding native units.
fprintf("\nProblem 6 no-step diagnostics:\n");
fprintf("P4 abstract state: dynamic/static estimate difference RMSE %.4g; gain gap %.2e.\n", sqrt(mean((resultP4NoStep.dynamic.xf - resultP4NoStep.static.xf).^2, "all")), gainRelativeErrorP4);
fprintf("P4 output vs physical truth: dynamic %.4g cm, static %.4g cm.\n", sqrt(mean((modelP4.C * resultP4NoStep.dynamic.xf - trueOutputNoStep6).^2, "all")), sqrt(mean((modelP4.C * resultP4NoStep.static.xf - trueOutputNoStep6).^2, "all")));
p5NoStepTankDynamicRMSE6 = sqrt(mean((tankHeightsFromMasses(resultP5NoStep.dynamic.xf(1:4, :), p6) - trueTankNoStep6).^2, "all"));
p5NoStepTankStaticRMSE6 = sqrt(mean((tankHeightsFromMasses(resultP5NoStep.static.xf(1:4, :), p6) - trueTankNoStep6).^2, "all"));
p5NoStepDisturbanceDynamicRMSE6 = sqrt(mean((resultP5NoStep.dynamic.xf(5:6, :) - trueDisturbanceNoStep6).^2, "all"));
p5NoStepDisturbanceStaticRMSE6 = sqrt(mean((resultP5NoStep.static.xf(5:6, :) - trueDisturbanceNoStep6).^2, "all"));
fprintf("P5 tank-height RMSE dynamic/static %.4g/%.4g cm; disturbance RMSE %.4g/%.4g cm^3/s.\n", p5NoStepTankDynamicRMSE6, p5NoStepTankStaticRMSE6, p5NoStepDisturbanceDynamicRMSE6, p5NoStepDisturbanceStaticRMSE6);
fprintf("P5 output vs physical truth: dynamic %.4g cm, static %.4g cm; gain gap %.2e.\n", sqrt(mean((modelP5.C * resultP5NoStep.dynamic.xf - trueOutputNoStep6).^2, "all")), sqrt(mean((modelP5.C * resultP5NoStep.static.xf - trueOutputNoStep6).^2, "all")), gainRelativeErrorP5);

% Plot P4 abstract states without interpreting them as physical tank heights.
plotP6EstimatePair(sampleTimes6, resultP4NoStep.static, resultP4NoStep.dynamic, 1:4, ["Abstract 1", "Abstract 2", "Abstract 3", "Abstract 4"], "P4 no-step: abstract estimates", fullfile(projectRoot, "figures", "problem6_6_p4_no_step_estimates.png"));
% Convert physical tank mass deviations to height deviations for tank-only figures.
p5NoStepStaticPlot6 = resultP5NoStep.static;
p5NoStepDynamicPlot6 = resultP5NoStep.dynamic;
p5NoStepStaticPlot6.xf(1:4, :) = tankHeightsFromMasses(p5NoStepStaticPlot6.xf(1:4, :), p6);
p5NoStepDynamicPlot6.xf(1:4, :) = tankHeightsFromMasses(p5NoStepDynamicPlot6.xf(1:4, :), p6);
tankOrder6 = [3, 4, 1, 2];
tankLabels6 = ["Tank 3 height dev. [cm]", "Tank 4 height dev. [cm]", "Tank 1 height dev. [cm]", "Tank 2 height dev. [cm]"];
plotP6TruthPair(sampleTimes6, trueTankNoStep6(tankOrder6, :), p5NoStepStaticPlot6.xf(tankOrder6, :), p5NoStepDynamicPlot6.xf(tankOrder6, :), tankLabels6, "P5 no-step: tank height dev.", fullfile(projectRoot, "figures", "problem6_6_p5_no_step_tank_truth.png"));
plotP6TruthPair(sampleTimes6, trueDisturbanceNoStep6, resultP5NoStep.static.xf(5:6, :), resultP5NoStep.dynamic.xf(5:6, :), ["D_3 [cm^3/s]", "D_4 [cm^3/s]"], "P5 no-step: disturbances", fullfile(projectRoot, "figures", "problem6_6_p5_no_step_disturbance_truth.png"));

% Compare dynamic gain convergence to the steady-state gains of both models.
figure("Visible", "off");
semilogy(sampleTimes6 / 60, max(gainGapP4, eps), "LineWidth", 1.2);
hold on;
semilogy(sampleTimes6 / 60, max(gainGapP5, eps), "LineWidth", 1.2);
grid on;
xlabel("Time [min]");
ylabel("||K_{dynamic} - K_{static}||_F");
legend("P4", "P5", "Location", "best");
title("No-step dynamic Kalman gains converge to static gains");
saveas(gcf, fullfile(projectRoot, "figures", "problem6_6_no_step_gain_convergence.png"));
close(gcf);

% Compare predicted outputs with physical truth and sampled measurements.
plotP6OutputFit(sampleTimes6, trueOutputNoStep6, measurementNoStep6, modelP4.C, resultP4NoStep.static, resultP4NoStep.dynamic, "P4 no-step output fit", fullfile(projectRoot, "figures", "problem6_6_p4_no_step_outputs.png"));
plotP6OutputFit(sampleTimes6, trueOutputNoStep6, measurementNoStep6, modelP5.C, resultP5NoStep.static, resultP5NoStep.dynamic, "P5 no-step output fit", fullfile(projectRoot, "figures", "problem6_6_p5_no_step_outputs.png"));

fprintf("Finished P6.6 Nonlinear model, no step\n");
%% 6.7 Offset-free filter models and step setup
% Apply the pump-disturbance mean change at the midpoint of the 30-minute run.
stepTime6 = 15 * 60;
stepMean6 = [80; 100];
meanSchedule6.times = [sim6.t0; stepTime6];
meanSchedule6.values = [p6.FMean, stepMean6];
rng(sim.randomSeed);

truthStepFine = simulateSDE_rerun(p6, sim6, meanSchedule6);
truthStep6 = truthStepFine.x(sampleIndex6, :);
truthStepD6 = truthStepFine.d(sampleIndex6, :);
measurementStep6 = truthStepFine.y(sampleIndex6, :).'- operpt6.y;
trueTankStep6 = tankHeightsFromMasses(truthStep6(:, 1:4).', p6) - operpt6.h;
trueDisturbanceStep6 = truthStepD6.' - p6.FMean;
trueOutputStep6 = truthStep6(:, 1:2).' ./ (p6.rho .* p6.A(1:2)) - operpt6.y;
trueMeanShift6 = zeros(2, sim6.N);
trueMeanShift6(:, sampleTimes6 >= stepTime6) = repmat(stepMean6 - p6.FMean, 1, sum(sampleTimes6 >= stepTime6));
stepSampleIndex6 = find(sampleTimes6 == stepTime6, 1);
assert(~isempty(stepSampleIndex6) && isequal(truthNoStepFine.v(sampleIndex6, :), truthStepFine.v(sampleIndex6, :)), "No-step and stepped experiments must use matching seeded measurement-noise samples.");

% Build the P4 augmentation with random-walk states for effective input biases.
% P4: integrating effective pump-input biases (not estimates of physical F3/F4).
modelP4N.A = [A4id, B4id; zeros(2, 4), eye(2)];
modelP4N.B = [B4id; zeros(2, 2)];
modelP4N.G = eye(6);
modelP4N.C = [C4id, zeros(2, 2)];
modelP4N.R = measurementCovariance6;
modelP4N.S = zeros(6, 2);
modelP4N.name = "P4 input-bias augmented";

% P5: continuous mean-shift augmentation from Chapter 6 equation (p6-model-p5-N).
Gamma6 = [zeros(4, 2); diag(p6.disturbanceRate)];
A8c = [aug6.A, Gamma6; zeros(2, 8)];
B8c = [aug6.B; zeros(2, 2)];

% Compare a baseline offset-noise scale with one decade below and above it.
SigmaNValues6 = [5];             % Per-sample standard deviations, can run comparison across many [cm^3/s]
p6Results.p4.linear.step = cell(1, numel(SigmaNValues6));
p6Results.p5.linear.step = cell(1, numel(SigmaNValues6));
p6Results.nonlinear.step = cell(1, numel(SigmaNValues6));
modelP5N = cell(1, numel(SigmaNValues6));
for noiseIndex = 1:numel(SigmaNValues6)
    sigmaN = SigmaNValues6(noiseIndex);
    SigmaNc = (sigmaN / sqrt(sim6.Ts)) * eye(2);
    G8c = [aug6.G, zeros(6, 2);
           zeros(2, 2), SigmaNc];
    [A8bar, B8bar, Q8bar] = c2dsde(A8c, B8c, G8c, sim6.Ts);
    Q8bar = (Q8bar + Q8bar') / 2;
    modelP5N{noiseIndex}.A = A8bar;
    modelP5N{noiseIndex}.B = B8bar;
    modelP5N{noiseIndex}.G = eye(8);
    modelP5N{noiseIndex}.C = [modelP5.C, zeros(2, 2)];
    modelP5N{noiseIndex}.Q = Q8bar;
    modelP5N{noiseIndex}.R = measurementCovariance6;
    modelP5N{noiseIndex}.S = zeros(8, 2);
    modelP5N{noiseIndex}.P0 = P0Scale6 * eye(8);
    modelP5N{noiseIndex}.name = "P5 mean-shift augmented";
    assertKalmanModelP6(modelP5N{noiseIndex});
    assert(norm(modelP5N{noiseIndex}.A(7:8, 7:8) - eye(2), "fro") < 1e-12, "P5 offset states must be discrete random walks.");
end

fprintf("Finished P6.7 Offset-free filter models and step setup\n");
%% 6.8 Problem 4 and 5 linear models, step
% Evaluate linear and nonlinear filters across the three offset-noise scales.
stepTransient = sampleTimes6 >= stepTime6 & sampleTimes6 < stepTime6 + 180;
stepSettled = sampleTimes6 >= stepTime6 + 180;
postStepSummary6 = zeros(numel(SigmaNValues6), 12);
linearStepSummary6 = zeros(numel(SigmaNValues6), 5);
for noiseIndex = 1:numel(SigmaNValues6)
    sigmaN = SigmaNValues6(noiseIndex);
    % Update the P4 augmented process covariance for the current sigma.
    modelP4N.Q = blkdiag(modelP4.Q, sigmaN^2 * eye(2));
    modelP4N.P0 = P0Scale6 * eye(6);
    assertKalmanModelP6(modelP4N);

    % The P4 step is an effective input bias, not a physical F3/F4 estimate.
    % A one-sample impulse sets the augmented bias state at the step.
    stepImpulse6 = zeros(2, sim6.N);
    stepImpulse6(:, stepSampleIndex6 - 1) = stepMean6 - p6.FMean;

    % Generate the P4 linear step truth and filter estimates using common noise.
    rng(sim.randomSeed + 620);
    modelP4NTruth = modelP4N;
    modelP4NTruth.P0 = zeros(6);
    modelP4NTruth.B = [modelP4N.B, [zeros(4, 2); eye(2)]];
    inputP4NTruth = [inputDeviation6; stepImpulse6];
    truthP4LinearStep = simulateLinearSystem(modelP4NTruth, inputP4NTruth, sim6.N);
    p4LinearStep.dynamic = kalmanFilter(modelP4N, truthP4LinearStep.y, inputDeviation6, "dynamic");
    p4LinearStep.static = kalmanFilter(modelP4N, truthP4LinearStep.y, inputDeviation6, "static");

    % Generate the corresponding P5 linear step realization for this sigma.
    rng(sim.randomSeed + 621);
    modelP5NTruth = modelP5N{noiseIndex};
    modelP5NTruth.P0 = zeros(8);
    modelP5NTruth.B = [modelP5NTruth.B, [zeros(6, 2); eye(2)]];
    inputP5NTruth = [inputDeviation6; stepImpulse6];
    truthP5LinearStep = simulateLinearSystem(modelP5NTruth, inputP5NTruth, sim6.N);
    p5LinearStep.dynamic = kalmanFilter(modelP5N{noiseIndex}, truthP5LinearStep.y, inputDeviation6, "dynamic");
    p5LinearStep.static = kalmanFilter(modelP5N{noiseIndex}, truthP5LinearStep.y, inputDeviation6, "static");

    % Retain linear-model truth, filter estimates, and sigma-specific scores.
    p6Results.p4.linear.step{noiseIndex}.truth = truthP4LinearStep;
    p6Results.p4.linear.step{noiseIndex}.filters = p4LinearStep;
    p6Results.p5.linear.step{noiseIndex}.truth = truthP5LinearStep;
    p6Results.p5.linear.step{noiseIndex}.filters = p5LinearStep;
    p6Results.p4.linear.step{noiseIndex}.sigmaN = sigmaN;
    p6Results.p5.linear.step{noiseIndex}.sigmaN = sigmaN;

    p4LinearOutput = modelP4N.C * truthP4LinearStep.x;
    p5LinearOutput = modelP5N{noiseIndex}.C * truthP5LinearStep.x;
    % Compute transient output RMSE for the final model comparison.
    p5LinearStepRMSE = [sqrt(mean((modelP5N{noiseIndex}.C * p5LinearStep.static.xp(:, stepTransient) - p5LinearOutput(:, stepTransient)).^2, "all")), sqrt(mean((modelP5N{noiseIndex}.C * p5LinearStep.dynamic.xp(:, stepTransient) - p5LinearOutput(:, stepTransient)).^2, "all"))];
    p4LinearStepRMSE = [sqrt(mean((modelP4N.C * p4LinearStep.static.xp(:, stepTransient) - p4LinearOutput(:, stepTransient)).^2, "all")), sqrt(mean((modelP4N.C * p4LinearStep.dynamic.xp(:, stepTransient) - p4LinearOutput(:, stepTransient)).^2, "all"))];
    linearStepSummary6(noiseIndex, :) = [sigmaN, p4LinearStepRMSE, p5LinearStepRMSE];
    p6Results.p4.linear.step{noiseIndex}.transientOutputRMSE = p4LinearStepRMSE;
    p6Results.p5.linear.step{noiseIndex}.transientOutputRMSE = p5LinearStepRMSE;

    sigmaTag = strrep(num2str(sigmaN), ".", "p");
    % Save linear-model output-fit figures under their section and sigma names.
    plotP6OutputFit(sampleTimes6, p4LinearOutput, truthP4LinearStep.y, modelP4N.C, p4LinearStep.static, p4LinearStep.dynamic, "P4 linear step output", fullfile(projectRoot, "figures", "problem6_8_p4_linear_step_outputs_sigma_" + sigmaTag + ".png"));
    plotP6OutputFit(sampleTimes6, p5LinearOutput, truthP5LinearStep.y, modelP5N{noiseIndex}.C, p5LinearStep.static, p5LinearStep.dynamic, "P5 linear step output", fullfile(projectRoot, "figures", "problem6_8_p5_linear_step_outputs_sigma_" + sigmaTag + ".png"));

    % Nonlinear model, step: evaluate both augmented filters on shared truth
    p4Step.dynamic = kalmanFilter(modelP4N, measurementStep6, inputDeviation6, "dynamic");
    p4Step.static = kalmanFilter(modelP4N, measurementStep6, inputDeviation6, "static");
    p5Step.dynamic = kalmanFilter(modelP5N{noiseIndex}, measurementStep6, inputDeviation6, "dynamic");
    p5Step.static = kalmanFilter(modelP5N{noiseIndex}, measurementStep6, inputDeviation6, "static");
    p6Results.nonlinear.step{noiseIndex}.truth = truthStepFine;
    p6Results.nonlinear.step{noiseIndex}.filters.p4 = p4Step;
    p6Results.nonlinear.step{noiseIndex}.filters.p5 = p5Step;
    p6Results.nonlinear.step{noiseIndex}.sigmaN = sigmaN;
    assert(all(isfinite(p4Step.dynamic.xf(:))) && all(isfinite(p5Step.dynamic.xf(:))), "Stepped filter returned nonfinite state estimates.");

    % Compare finite-horizon gains with each augmented model's steady-state gain.
    [~, staticGainP4N] = kalmanStatic(modelP4N);
    [~, staticGainP5N] = kalmanStatic(modelP5N{noiseIndex});
    % Zero dummy data isolate covariance-driven gain convergence.
    gainProbeSamples6 = 3000;
    gainProbeInput6 = zeros(2, gainProbeSamples6);
    gainProbeOutput6 = zeros(2, gainProbeSamples6);
    gainProbeP4N = kalmanFilter(modelP4N, gainProbeOutput6, gainProbeInput6, "dynamic");
    gainProbeP5N = kalmanFilter(modelP5N{noiseIndex}, gainProbeOutput6, gainProbeInput6, "dynamic");
    probeGapP4N = norm(gainProbeP4N.K(:, :, end) - staticGainP4N, "fro") / max(norm(staticGainP4N, "fro"), eps);
    probeGapP5N = norm(gainProbeP5N.K(:, :, end) - staticGainP5N, "fro") / max(norm(staticGainP5N, "fro"), eps);
    assert(probeGapP4N < 1e-5 && probeGapP5N < 1e-5, "Dynamic augmented-model gains did not converge to their static gains.");
    p4StepGainGap = norm(p4Step.dynamic.K(:, :, end) - staticGainP4N, "fro") / max(norm(staticGainP4N, "fro"), eps);
    p5StepGainGap = norm(p5Step.dynamic.K(:, :, end) - staticGainP5N, "fro") / max(norm(staticGainP5N, "fro"), eps);

    % Calculate output predictions and summary metrics for both augmented models.
    p4PredictionDynamic = modelP4N.C * p4Step.dynamic.xp;
    p4PredictionStatic = modelP4N.C * p4Step.static.xp;
    p5PredictionDynamic = modelP5N{noiseIndex}.C * p5Step.dynamic.xp;
    p5PredictionStatic = modelP5N{noiseIndex}.C * p5Step.static.xp;
    postStepSummary6(noiseIndex, 1) = sigmaN;
    postStepSummary6(noiseIndex, 2) = sqrt(mean((p5PredictionDynamic(:, stepTransient) - trueOutputStep6(:, stepTransient)).^2, "all"));
    postStepSummary6(noiseIndex, 3) = max(abs(mean(p5PredictionDynamic(:, stepSettled) - trueOutputStep6(:, stepSettled), 2)));
    postStepSummary6(noiseIndex, 4) = sqrt(mean((p5Step.dynamic.xf(7:8, stepSettled) - trueMeanShift6(:, stepSettled)).^2, "all"));
    postStepSummary6(noiseIndex, 5) = sqrt(mean((p5PredictionStatic(:, stepTransient) - trueOutputStep6(:, stepTransient)).^2, "all"));
    postStepSummary6(noiseIndex, 6) = max(abs(mean(p5PredictionStatic(:, stepSettled) - trueOutputStep6(:, stepSettled), 2)));
    postStepSummary6(noiseIndex, 7) = sqrt(mean((p5Step.static.xf(7:8, stepSettled) - trueMeanShift6(:, stepSettled)).^2, "all"));
    postStepSummary6(noiseIndex, 8) = sqrt(mean((p4PredictionDynamic(:, stepTransient) - trueOutputStep6(:, stepTransient)).^2, "all"));
    postStepSummary6(noiseIndex, 9) = max(abs(mean(p4PredictionDynamic(:, stepSettled) - trueOutputStep6(:, stepSettled), 2)));
    postStepSummary6(noiseIndex, 10) = sqrt(mean((p4PredictionStatic(:, stepTransient) - trueOutputStep6(:, stepTransient)).^2, "all"));
    postStepSummary6(noiseIndex, 11) = max(abs(mean(p4PredictionStatic(:, stepSettled) - trueOutputStep6(:, stepSettled), 2)));
    postStepSummary6(noiseIndex, 12) = max(p4StepGainGap, p5StepGainGap);

    % Confirm completion of each sigma case; detailed metrics are summarized below.
    fprintf("Completed nonlinear step estimates for sigma_N=%g cm^3/s.\n", sigmaN);

    % Convert physical tank states to height deviations only in plot copies.
    p5StepStaticPlot6 = p5Step.static;
    p5StepDynamicPlot6 = p5Step.dynamic;
    p5StepStaticPlot6.xf(1:4, :) = tankHeightsFromMasses(p5StepStaticPlot6.xf(1:4, :), p6);
    p5StepDynamicPlot6.xf(1:4, :) = tankHeightsFromMasses(p5StepDynamicPlot6.xf(1:4, :), p6);
    % Plot the P4 augmented state estimates and output fit.
    plotP6EstimatePair(sampleTimes6, p4Step.static, p4Step.dynamic, 1:6, ["Abstract 1", "Abstract 2", "Abstract 3", "Abstract 4", "Eff. input bias 1", "Eff. input bias 2"], "P4 step: estimates", fullfile(projectRoot, "figures", "problem6_8_p4_step_estimates_sigma_" + sigmaTag + ".png"));
    plotP6OutputFit(sampleTimes6, trueOutputStep6, measurementStep6, modelP4N.C, p4Step.static, p4Step.dynamic, "P4 step output", fullfile(projectRoot, "figures", "problem6_8_p4_step_outputs_sigma_" + sigmaTag + ".png"));

    % Compare tank-height and D/N truth with both filters in separate figures.
    plotP6TruthPair(sampleTimes6, trueTankStep6(tankOrder6, :), p5StepStaticPlot6.xf(tankOrder6, :), p5StepDynamicPlot6.xf(tankOrder6, :), tankLabels6, "P5 step: tank height dev.", fullfile(projectRoot, "figures", "problem6_8_p5_step_tank_heights_sigma_" + sigmaTag + ".png"));
    plotP6TruthPair(sampleTimes6, [trueDisturbanceStep6; trueMeanShift6], p5Step.static.xf(5:8, :), p5Step.dynamic.xf(5:8, :), ["D_3 [cm^3/s]", "D_4 [cm^3/s]", "N_3 [cm^3/s]", "N_4 [cm^3/s]"], "P5 step: D and N estimates", fullfile(projectRoot, "figures", "problem6_8_p5_step_D_N_sigma_" + sigmaTag + ".png"));
    plotP6OutputFit(sampleTimes6, trueOutputStep6, measurementStep6, modelP5N{noiseIndex}.C, p5Step.static, p5Step.dynamic, "P5 step output", fullfile(projectRoot, "figures", "problem6_8_p5_step_outputs_sigma_" + sigmaTag + ".png"));
end

% Report settled output bias, P5 offset-state error, and the gain-convergence gap.
fprintf("\nPost-step summary (settled interval begins 180 s after the step):\n");
postStepSummaryTable6 = array2table(postStepSummary6(:, [1, 3, 6, 9, 11, 4, 7, 12]), "VariableNames", ["sigmaN_cm3s", "P5_dynamic_settled_bias_cm", "P5_static_settled_bias_cm", "P4_dynamic_settled_bias_cm", "P4_static_settled_bias_cm", "P5_dynamic_offset_RMSE_cm3s", "P5_static_offset_RMSE_cm3s", "max_dynamic_static_gain_gap"]);
disp(postStepSummaryTable6);
linearStepSummaryTable6 = array2table(linearStepSummary6, "VariableNames", ["sigmaN_cm3s", "P4_static_RMSE_cm", "P4_dynamic_RMSE_cm", "P5_static_RMSE_cm", "P5_dynamic_RMSE_cm"]);
fprintf("\nLinear-model step transient output RMSE [cm] by sigma_N:\n");
fprintf("Each row is one sigma_N case; RMSE compares each filter with its matching linear-model truth over the 180-second post-step transient.\n");
disp(linearStepSummaryTable6);

% Mark completion of the linear and nonlinear stepped-filter comparisons.
fprintf("Finished P6.8 Problem 4 and 5 linear models, step\n");
%% 6.9 Compare and plot output estimation errors
% Compare no-step and baseline-sigma step output errors across both model types.
% Aggregate no-step output errors for linear and nonlinear truth experiments.
nonlinearNoStepRMSE6 = [sqrt(mean((modelP4.C * resultP4NoStep.static.xf - trueOutputNoStep6).^2, "all")), sqrt(mean((modelP4.C * resultP4NoStep.dynamic.xf - trueOutputNoStep6).^2, "all")); sqrt(mean((modelP5.C * resultP5NoStep.static.xf - trueOutputNoStep6).^2, "all")), sqrt(mean((modelP5.C * resultP5NoStep.dynamic.xf - trueOutputNoStep6).^2, "all"))];
noStepRMSEComparison6 = [p6Results.p4.linear.noStep.outputRMSE; p6Results.p5.linear.noStep.outputRMSE; nonlinearNoStepRMSE6];

% Select the baseline for the nonlinear step comparison.
middleNoiseIndex6 = find(SigmaNValues6 == 5, 1);
% Compare transient step output RMSE at the baseline sigma across model types.
stepRMSEComparison6 = [linearStepSummary6(middleNoiseIndex6, 2:3); linearStepSummary6(middleNoiseIndex6, 4:5); postStepSummary6(middleNoiseIndex6, [10, 8]); postStepSummary6(middleNoiseIndex6, [5, 2])];
% Render both RMSE comparisons with a shared filter legend.
figure("Visible", "off");
tiledlayout(1, 2, "TileSpacing", "compact");
nexttile;
bar(noStepRMSEComparison6);
grid on;
set(gca, "XTickLabel", ["P4 linear", "P5 linear", "P4 on nonlinear", "P5 on nonlinear"]);
ylabel("Output RMSE [cm]");
legend("Static", "Dynamic", "Location", "best");
title("No-step output estimation");
nexttile;
bar(stepRMSEComparison6);
grid on;
set(gca, "XTickLabel", ["P4 linear", "P5 linear", "P4 on nonlinear", "P5 on nonlinear"]);
ylabel("Transient output RMSE [cm]");
title("Step output estimation, \sigma_N = 2");
exportgraphics(gcf, fullfile(projectRoot, "figures", "problem6_9_output_rmse_comparison.png"), "Resolution", 400);
close(gcf);

fprintf("Finished P6.9 Compare and plot output estimation errors\n");