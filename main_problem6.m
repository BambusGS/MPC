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
sim6.tf = 45 * 60;
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
tankOrder6 = [3, 4, 1, 2];
tankLabels6 = ["Tank 3 height dev. [cm]", "Tank 4 height dev. [cm]", "Tank 1 height dev. [cm]", "Tank 2 height dev. [cm]"];

fprintf("Finished P6.1 Setup\n");
%% 6.2 P4 10% identified model and four-state realization

% Load CH4 markov parameters
loadedMarkov = load(fullfile(projectRoot, "markov.mat"), "markov");
H4transfer = loadedMarkov.markov(:, :, :, 1); % 10% step-size model

hankelSizeP4 = 44;
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

% Perform P5 model linearization and augmentation.
lin6 = computeJacobians(operpt6, p6);
aug6 = computeAugmentedModel(lin6, p6);

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

%% 6.4 No-step: P4 linear-model truth vs. static/dynamic P4 Kalman filters
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
plotP6OutputFit(sampleTimes6, modelP4.C * truthP4LinearNoStep.x, truthP4LinearNoStep.y, modelP4.C, resultP4LinearNoStep.static, resultP4LinearNoStep.dynamic, "P4 original linear no-step", "Linear truth", fullfile(projectRoot, "figures", "problem6_4_p4_linear_no_step_outputs.png"));
plotP6TruthPair(sampleTimes6, truthP4LinearNoStep.x, resultP4LinearNoStep.static.xf, resultP4LinearNoStep.dynamic.xf, ["Abstract 1", "Abstract 2", "Abstract 3", "Abstract 4"], "P4 original no-step: abstract states", "Linear truth", fullfile(projectRoot, "figures", "problem6_4_p4_linear_no_step_abstract_states.png"));

fprintf("Finished P6.4 P4 linear model, no step\n");
%% 6.5 No-step: P5 linear-model truth vs. static/dynamic P5 Kalman filters
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
plotP6OutputFit(sampleTimes6, modelP5.C * truthP5LinearNoStep.x, truthP5LinearNoStep.y, modelP5.C, resultP5LinearNoStep.static, resultP5LinearNoStep.dynamic, "P5 original linear no-step", "Linear truth", fullfile(projectRoot, "figures", "problem6_5_p5_linear_no_step_outputs.png"));
p5LinearTruthTankHeight6 = tankHeightsFromMasses(truthP5LinearNoStep.x(1:4, :), p6);
p5LinearStaticTankHeight6 = tankHeightsFromMasses(resultP5LinearNoStep.static.xf(1:4, :), p6);
p5LinearDynamicTankHeight6 = tankHeightsFromMasses(resultP5LinearNoStep.dynamic.xf(1:4, :), p6);
plotP6TruthPair(sampleTimes6, p5LinearTruthTankHeight6(tankOrder6, :), p5LinearStaticTankHeight6(tankOrder6, :), p5LinearDynamicTankHeight6(tankOrder6, :), tankLabels6, "P5 original no-step: tank height dev.", "Linear truth", fullfile(projectRoot, "figures", "problem6_5_p5_linear_no_step_tank_heights.png"));
plotP6TruthPair(sampleTimes6, truthP5LinearNoStep.x(5:6, :), resultP5LinearNoStep.static.xf(5:6, :), resultP5LinearNoStep.dynamic.xf(5:6, :), ["D_3 [cm^3/s]", "D_4 [cm^3/s]"], "P5 original no-step: disturbances", "Linear truth", fullfile(projectRoot, "figures", "problem6_5_p5_linear_no_step_disturbances.png"));

fprintf("Finished P6.5 P5 linear model, no step\n");
%% 6.6 No-step: nonlinear SDE truth vs. static/dynamic P4/P5 linear Kalman filters
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

% Convert physical tank mass deviations to height deviations for tank-only figures.
p5NoStepStaticPlot6 = resultP5NoStep.static;
p5NoStepDynamicPlot6 = resultP5NoStep.dynamic;
p5NoStepStaticPlot6.xf(1:4, :) = tankHeightsFromMasses(p5NoStepStaticPlot6.xf(1:4, :), p6);
p5NoStepDynamicPlot6.xf(1:4, :) = tankHeightsFromMasses(p5NoStepDynamicPlot6.xf(1:4, :), p6);
plotP6TruthPair(sampleTimes6, trueTankNoStep6(tankOrder6, :), p5NoStepStaticPlot6.xf(tankOrder6, :), p5NoStepDynamicPlot6.xf(tankOrder6, :), tankLabels6, "P5 original no-step: tank height dev.", "Nonlinear truth", fullfile(projectRoot, "figures", "problem6_6_p5_no_step_tank_truth.png"));
plotP6TruthPair(sampleTimes6, trueDisturbanceNoStep6, resultP5NoStep.static.xf(5:6, :), resultP5NoStep.dynamic.xf(5:6, :), ["D_3 [cm^3/s]", "D_4 [cm^3/s]"], "P5 original no-step: disturbances", "Nonlinear truth", fullfile(projectRoot, "figures", "problem6_6_p5_no_step_disturbance_truth.png"));

% Compare dynamic gain convergence to the steady-state gains of both models.
figure("Visible", "off");
semilogy(sampleTimes6 / 60, max(gainGapP4 / max(norm(staticGainP4, "fro"), eps), eps), "LineWidth", 1.2);
hold on;
semilogy(sampleTimes6 / 60, max(gainGapP5 / max(norm(staticGainP5, "fro"), eps), eps), "LineWidth", 1.2);
grid on;
xlabel("Time [min]");
ylabel("Relative ||K_{dynamic} - K_{static}||_F");
legend("P4 original", "P5 original", "Location", "best");
title("No-step dynamic Kalman gains converge to static gains");
saveas(gcf, fullfile(projectRoot, "figures", "problem6_6_no_step_gain_convergence.png"));
close(gcf);

% Compare predicted outputs with physical truth and sampled measurements.
plotP6OutputFit(sampleTimes6, trueOutputNoStep6, measurementNoStep6, modelP4.C, resultP4NoStep.static, resultP4NoStep.dynamic, "P4 original no-step output fit", "Nonlinear truth", fullfile(projectRoot, "figures", "problem6_6_p4_no_step_outputs.png"));
plotP6OutputFit(sampleTimes6, trueOutputNoStep6, measurementNoStep6, modelP5.C, resultP5NoStep.static, resultP5NoStep.dynamic, "P5 original no-step output fit", "Nonlinear truth", fullfile(projectRoot, "figures", "problem6_6_p5_no_step_outputs.png"));

fprintf("Finished P6.6 Nonlinear model, no step\n");
%% 6.7 Step scenario setup: nonlinear SDE truth and augmented P4/P5 linear filters
% Change the F3 mean at 15 min and the F4 mean at 30 min.
stepTimes6 = [15, 30] * 60;
stepDeltaValues6 = [30, 0;
                    0, -20];
stepMeanValues6 = [p6.FMean, p6.FMean + cumsum(stepDeltaValues6, 2)];
meanSchedule6.times = [sim6.t0; stepTimes6(:)];
meanSchedule6.values = stepMeanValues6;
rng(sim.randomSeed);

truthStepFine = simulateSDE_rerun(p6, sim6, meanSchedule6);
truthStep6 = truthStepFine.x(sampleIndex6, :);
truthStepD6 = truthStepFine.d(sampleIndex6, :);
measurementStep6 = truthStepFine.y(sampleIndex6, :).'- operpt6.y;
trueTankStep6 = tankHeightsFromMasses(truthStep6(:, 1:4).', p6) - operpt6.h;
trueDisturbanceStep6 = truthStepD6.' - p6.FMean;
trueOutputStep6 = truthStep6(:, 1:2).' ./ (p6.rho .* p6.A(1:2)) - operpt6.y;
trueMeanShift6 = zeros(2, sim6.N);
for stepIndex6 = 1:numel(stepTimes6)
    stepMask6 = sampleTimes6.' >= stepTimes6(stepIndex6);
    currentMeanShift6 = stepMeanValues6(:, stepIndex6 + 1) - p6.FMean;
    trueMeanShift6(:, stepMask6) = repmat(currentMeanShift6, 1, sum(stepMask6));
end
stepSampleIndices6 = zeros(1, numel(stepTimes6));
for stepIndex6 = 1:numel(stepTimes6)
    stepSampleIndex6 = find(sampleTimes6 == stepTimes6(stepIndex6), 1);
    assert(~isempty(stepSampleIndex6), "Every disturbance step must align with a filter sample.");
    stepSampleIndices6(stepIndex6) = stepSampleIndex6;
end
assert(isequal(truthNoStepFine.v(sampleIndex6, :), truthStepFine.v(sampleIndex6, :)), "No-step and stepped experiments must use matching seeded measurement-noise samples.");

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

% Each configured offset-noise scale generates its complete stepped-figure set.
SigmaNValues6 = [5];             % Per-sample standard deviations, can run comparison across many if needed [cm^3/s]
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
end

fprintf("Finished P6.7 Offset-free filter models and step setup\n");
%% 6.8 Step: linear-model truth and nonlinear SDE truth vs. augmented P4/P5 filters
% Evaluate linear and nonlinear filters at every configured offset-noise scale.
stepEventCount6 = numel(stepTimes6);
stepEventEndTimes6 = [stepTimes6(2:end), sim6.tf + sim6.Ts];
stepTransientMasks6 = false(stepEventCount6, sim6.N);
stepSettledMasks6 = false(stepEventCount6, sim6.N);
for stepIndex6 = 1:stepEventCount6
    stepTransientEnd6 = min(stepTimes6(stepIndex6) + 180, stepEventEndTimes6(stepIndex6));
    stepTransientMasks6(stepIndex6, :) = sampleTimes6.' >= stepTimes6(stepIndex6) & sampleTimes6.' < stepTransientEnd6;
    stepSettledMasks6(stepIndex6, :) = sampleTimes6.' >= stepTimes6(stepIndex6) + 180 & sampleTimes6.' < stepEventEndTimes6(stepIndex6);
end
stepSummaryRowCount6 = numel(SigmaNValues6) * stepEventCount6;
postStepSummary6 = zeros(stepSummaryRowCount6, 13);
linearStepSummary6 = zeros(stepSummaryRowCount6, 6);
for noiseIndex = 1:numel(SigmaNValues6)
    sigmaN = SigmaNValues6(noiseIndex);
    % Update the P4 augmented process covariance for the current sigma.
    modelP4N.Q = blkdiag(modelP4.Q, sigmaN^2 * eye(2));
    modelP4N.P0 = P0Scale6 * eye(6);
    assertKalmanModelP6(modelP4N);

    % P4 estimates effective input biases, not physical F3/F4 flows.
    % One-sample impulses set the augmented bias states at each step.
    stepImpulse6 = zeros(2, sim6.N);
    for stepIndex6 = 1:stepEventCount6
        stepImpulse6(:, stepSampleIndices6(stepIndex6) - 1) = stepDeltaValues6(:, stepIndex6);
    end

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
    p4LinearStepRMSEByEvent6 = zeros(stepEventCount6, 2);
    p5LinearStepRMSEByEvent6 = zeros(stepEventCount6, 2);
    for stepIndex6 = 1:stepEventCount6
        stepTransientMask6 = stepTransientMasks6(stepIndex6, :);
        p5LinearStepRMSEByEvent6(stepIndex6, :) = [sqrt(mean((modelP5N{noiseIndex}.C * p5LinearStep.static.xp(:, stepTransientMask6) - p5LinearOutput(:, stepTransientMask6)).^2, "all")), sqrt(mean((modelP5N{noiseIndex}.C * p5LinearStep.dynamic.xp(:, stepTransientMask6) - p5LinearOutput(:, stepTransientMask6)).^2, "all"))];
        p4LinearStepRMSEByEvent6(stepIndex6, :) = [sqrt(mean((modelP4N.C * p4LinearStep.static.xp(:, stepTransientMask6) - p4LinearOutput(:, stepTransientMask6)).^2, "all")), sqrt(mean((modelP4N.C * p4LinearStep.dynamic.xp(:, stepTransientMask6) - p4LinearOutput(:, stepTransientMask6)).^2, "all"))];
        summaryRow6 = (noiseIndex - 1) * stepEventCount6 + stepIndex6;
        linearStepSummary6(summaryRow6, :) = [sigmaN, stepTimes6(stepIndex6) / 60, p4LinearStepRMSEByEvent6(stepIndex6, :), p5LinearStepRMSEByEvent6(stepIndex6, :)];
    end
    allStepTransientMask6 = any(stepTransientMasks6, 1);
    p4LinearStepRMSE = [sqrt(mean((modelP4N.C * p4LinearStep.static.xp(:, allStepTransientMask6) - p4LinearOutput(:, allStepTransientMask6)).^2, "all")), sqrt(mean((modelP4N.C * p4LinearStep.dynamic.xp(:, allStepTransientMask6) - p4LinearOutput(:, allStepTransientMask6)).^2, "all"))];
    p5LinearStepRMSE = [sqrt(mean((modelP5N{noiseIndex}.C * p5LinearStep.static.xp(:, allStepTransientMask6) - p5LinearOutput(:, allStepTransientMask6)).^2, "all")), sqrt(mean((modelP5N{noiseIndex}.C * p5LinearStep.dynamic.xp(:, allStepTransientMask6) - p5LinearOutput(:, allStepTransientMask6)).^2, "all"))];
    p6Results.p4.linear.step{noiseIndex}.transientOutputRMSE = p4LinearStepRMSE;
    p6Results.p4.linear.step{noiseIndex}.transientOutputRMSEByEvent = p4LinearStepRMSEByEvent6;
    p6Results.p5.linear.step{noiseIndex}.transientOutputRMSE = p5LinearStepRMSE;
    p6Results.p5.linear.step{noiseIndex}.transientOutputRMSEByEvent = p5LinearStepRMSEByEvent6;
    p6Results.p4.linear.step{noiseIndex}.stepTimes = stepTimes6;
    p6Results.p5.linear.step{noiseIndex}.stepTimes = stepTimes6;

    sigmaTag = strrep(num2str(sigmaN), ".", "p");
    % Save linear-model output-fit figures under their section and sigma names.
    plotP6OutputFit(sampleTimes6, p4LinearOutput, truthP4LinearStep.y, modelP4N.C, p4LinearStep.static, p4LinearStep.dynamic, "P4 augmented step output", "Linear truth", fullfile(projectRoot, "figures", "problem6_8_p4_augmented_linear_step_outputs_sigma_" + sigmaTag + ".png"), stepTimes6);
    plotP6OutputFit(sampleTimes6, p5LinearOutput, truthP5LinearStep.y, modelP5N{noiseIndex}.C, p5LinearStep.static, p5LinearStep.dynamic, "P5 augmented step output", "Linear truth", fullfile(projectRoot, "figures", "problem6_8_p5_augmented_linear_step_outputs_sigma_" + sigmaTag + ".png"), stepTimes6);
    plotP6TruthPair(sampleTimes6, truthP4LinearStep.x, p4LinearStep.static.xf, p4LinearStep.dynamic.xf, ["Abstract 1", "Abstract 2", "Abstract 3", "Abstract 4", "Dist. integrator 1", "Dist. integrator 2"], "P4 augmented step: abstract states", "Linear truth", fullfile(projectRoot, "figures", "problem6_8_p4_augmented_linear_step_abstract_states_sigma_" + sigmaTag + ".png"), stepTimes6);

    p5LinearStepStaticPlot6 = p5LinearStep.static;
    p5LinearStepDynamicPlot6 = p5LinearStep.dynamic;
    p5LinearStepStaticPlot6.xf(1:4, :) = tankHeightsFromMasses(p5LinearStepStaticPlot6.xf(1:4, :), p6);
    p5LinearStepDynamicPlot6.xf(1:4, :) = tankHeightsFromMasses(p5LinearStepDynamicPlot6.xf(1:4, :), p6);
    p5LinearStepTruthTankHeight6 = tankHeightsFromMasses(truthP5LinearStep.x(1:4, :), p6);
    plotP6TruthPair(sampleTimes6, p5LinearStepTruthTankHeight6(tankOrder6, :), p5LinearStepStaticPlot6.xf(tankOrder6, :), p5LinearStepDynamicPlot6.xf(tankOrder6, :), tankLabels6, "P5 augmented step: tank height dev.", "Linear truth", fullfile(projectRoot, "figures", "problem6_8_p5_augmented_linear_step_tank_heights_sigma_" + sigmaTag + ".png"), stepTimes6);
    plotP6TruthPair(sampleTimes6, truthP5LinearStep.x(5:8, :), p5LinearStep.static.xf(5:8, :), p5LinearStep.dynamic.xf(5:8, :), ["D_3 [cm^3/s]", "D_4 [cm^3/s]", "N_3 [cm^3/s]", "N_4 [cm^3/s]"], "P5 augmented step: disturbances and N estimates", "Linear truth", fullfile(projectRoot, "figures", "problem6_8_p5_augmented_linear_step_D_N_sigma_" + sigmaTag + ".png"), stepTimes6);

    % Nonlinear model, step: evaluate both augmented filters on shared truth
    p4Step.dynamic = kalmanFilter(modelP4N, measurementStep6, inputDeviation6, "dynamic");
    p4Step.static = kalmanFilter(modelP4N, measurementStep6, inputDeviation6, "static");
    p5Step.dynamic = kalmanFilter(modelP5N{noiseIndex}, measurementStep6, inputDeviation6, "dynamic");
    p5Step.static = kalmanFilter(modelP5N{noiseIndex}, measurementStep6, inputDeviation6, "static");
    p6Results.nonlinear.step{noiseIndex}.truth = truthStepFine;
    p6Results.nonlinear.step{noiseIndex}.filters.p4 = p4Step;
    p6Results.nonlinear.step{noiseIndex}.filters.p5 = p5Step;
    p6Results.nonlinear.step{noiseIndex}.sigmaN = sigmaN;

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
    gainGapP4Step6 = zeros(1, sim6.N);
    gainGapP5Step6 = zeros(1, sim6.N);
    for k = 1:sim6.N
        gainGapP4Step6(k) = norm(p4Step.dynamic.K(:, :, k) - staticGainP4N, "fro") / max(norm(staticGainP4N, "fro"), eps);
        gainGapP5Step6(k) = norm(p5Step.dynamic.K(:, :, k) - staticGainP5N, "fro") / max(norm(staticGainP5N, "fro"), eps);
    end

    % Calculate output predictions and summary metrics for both augmented models.
    p4PredictionDynamic = modelP4N.C * p4Step.dynamic.xp;
    p4PredictionStatic = modelP4N.C * p4Step.static.xp;
    p5PredictionDynamic = modelP5N{noiseIndex}.C * p5Step.dynamic.xp;
    p5PredictionStatic = modelP5N{noiseIndex}.C * p5Step.static.xp;
    for stepIndex6 = 1:stepEventCount6
        stepTransientMask6 = stepTransientMasks6(stepIndex6, :);
        stepSettledMask6 = stepSettledMasks6(stepIndex6, :);
        summaryRow6 = (noiseIndex - 1) * stepEventCount6 + stepIndex6;
        postStepSummary6(summaryRow6, :) = [sigmaN, stepTimes6(stepIndex6) / 60, ...
            sqrt(mean((p5PredictionDynamic(:, stepTransientMask6) - trueOutputStep6(:, stepTransientMask6)).^2, "all")), ...
            max(abs(mean(p5PredictionDynamic(:, stepSettledMask6) - trueOutputStep6(:, stepSettledMask6), 2))), ...
            sqrt(mean((p5Step.dynamic.xf(7:8, stepSettledMask6) - trueMeanShift6(:, stepSettledMask6)).^2, "all")), ...
            sqrt(mean((p5PredictionStatic(:, stepTransientMask6) - trueOutputStep6(:, stepTransientMask6)).^2, "all")), ...
            max(abs(mean(p5PredictionStatic(:, stepSettledMask6) - trueOutputStep6(:, stepSettledMask6), 2))), ...
            sqrt(mean((p5Step.static.xf(7:8, stepSettledMask6) - trueMeanShift6(:, stepSettledMask6)).^2, "all")), ...
            sqrt(mean((p4PredictionDynamic(:, stepTransientMask6) - trueOutputStep6(:, stepTransientMask6)).^2, "all")), ...
            max(abs(mean(p4PredictionDynamic(:, stepSettledMask6) - trueOutputStep6(:, stepSettledMask6), 2))), ...
            sqrt(mean((p4PredictionStatic(:, stepTransientMask6) - trueOutputStep6(:, stepTransientMask6)).^2, "all")), ...
            max(abs(mean(p4PredictionStatic(:, stepSettledMask6) - trueOutputStep6(:, stepSettledMask6), 2))), ...
            max(p4StepGainGap, p5StepGainGap)];
    end

    % Confirm completion of each sigma case; detailed metrics are summarized below.
    fprintf("Completed nonlinear step estimates for sigma_N=%g cm^3/s.\n", sigmaN);

    % Convert physical tank states to height deviations only in plot copies.
    p5StepStaticPlot6 = p5Step.static;
    p5StepDynamicPlot6 = p5Step.dynamic;
    p5StepStaticPlot6.xf(1:4, :) = tankHeightsFromMasses(p5StepStaticPlot6.xf(1:4, :), p6);
    p5StepDynamicPlot6.xf(1:4, :) = tankHeightsFromMasses(p5StepDynamicPlot6.xf(1:4, :), p6);
    plotP6OutputFit(sampleTimes6, trueOutputStep6, measurementStep6, modelP4N.C, p4Step.static, p4Step.dynamic, "P4 augmented step output", "Nonlinear truth", fullfile(projectRoot, "figures", "problem6_8_p4_augmented_nonlinear_step_outputs_sigma_" + sigmaTag + ".png"), stepTimes6);

    % Compare tank-height and D/N truth with both filters in separate figures.
    plotP6TruthPair(sampleTimes6, trueTankStep6(tankOrder6, :), p5StepStaticPlot6.xf(tankOrder6, :), p5StepDynamicPlot6.xf(tankOrder6, :), tankLabels6, "P5 augmented step: tank height dev.", "Nonlinear truth", fullfile(projectRoot, "figures", "problem6_8_p5_augmented_nonlinear_step_tank_heights_sigma_" + sigmaTag + ".png"), stepTimes6);
    plotP6TruthPair(sampleTimes6, [trueDisturbanceStep6; trueMeanShift6], p5Step.static.xf(5:8, :), p5Step.dynamic.xf(5:8, :), ["D_3 [cm^3/s]", "D_4 [cm^3/s]", "N_3 [cm^3/s]", "N_4 [cm^3/s]"], "P5 augmented step: disturbances and N estimates", "Nonlinear truth", fullfile(projectRoot, "figures", "problem6_8_p5_augmented_nonlinear_step_D_N_sigma_" + sigmaTag + ".png"), stepTimes6);
    plotP6OutputFit(sampleTimes6, trueOutputStep6, measurementStep6, modelP5N{noiseIndex}.C, p5Step.static, p5Step.dynamic, "P5 augmented step output", "Nonlinear truth", fullfile(projectRoot, "figures", "problem6_8_p5_augmented_nonlinear_step_outputs_sigma_" + sigmaTag + ".png"), stepTimes6);

    figure("Visible", "off");
    semilogy(sampleTimes6 / 60, max(gainGapP4Step6, eps), "LineWidth", 1.2);
    hold on;
    semilogy(sampleTimes6 / 60, max(gainGapP5Step6, eps), "LineWidth", 1.2);
    grid on;
    xlabel("Time [min]");
    ylabel("Relative ||K_{dynamic} - K_{static}||_F");
    legend("P4 augmented", "P5 augmented", "Location", "best");
    title("Step dynamic Kalman gains converge to static gains");
    saveas(gcf, fullfile(projectRoot, "figures", "problem6_8_augmented_gain_convergence_sigma_" + sigmaTag + ".png"));
    close(gcf);
end

% Report transient and settled metrics separately for each step event.
fprintf("\nPost-step summary by event (settled interval starts 180 s after each step):\n");
postStepSummaryTable6 = array2table(postStepSummary6, "VariableNames", ["sigmaN_cm3s", "stepTime_min", "P5_dynamic_transient_RMSE_cm", "P5_dynamic_settled_bias_cm", "P5_dynamic_offset_RMSE_cm3s", "P5_static_transient_RMSE_cm", "P5_static_settled_bias_cm", "P5_static_offset_RMSE_cm3s", "P4_dynamic_transient_RMSE_cm", "P4_dynamic_settled_bias_cm", "P4_static_transient_RMSE_cm", "P4_static_settled_bias_cm", "max_dynamic_static_gain_gap"]);
disp(postStepSummaryTable6);
linearStepSummaryTable6 = array2table(linearStepSummary6, "VariableNames", ["sigmaN_cm3s", "stepTime_min", "P4_static_RMSE_cm", "P4_dynamic_RMSE_cm", "P5_static_RMSE_cm", "P5_dynamic_RMSE_cm"]);
fprintf("\nLinear-model step transient output RMSE [cm] by event and sigma_N (first 180 s after each step):\n");
fprintf("Each row compares each filter with its matching model-specific linear truth during that event's transient window.\n");
disp(linearStepSummaryTable6);

% Mark completion of the linear and nonlinear stepped-filter comparisons.
fprintf("Finished P6.8 Problem 4 and 5 linear models, step\n");