%% Problem 6 setup
projectRoot = fileparts(mfilename("fullpath"));
addpath(fullfile(projectRoot, "functions"));
addpath(fullfile(projectRoot, "Problem6_KF_Demo", "functions"));
p = createParameters();
sim = createSimSettings(p);

%% 6.1 Configuration, operating point, and repeatable simulations
p6 = p;
p6Results = struct();
p6.FMean = [50; 100];                     % Nominal [F3; F4] [cm^3/s]
p6.disturbanceRate = [0.05; 0.05];         % Mean-reversion rates [1/s]
p6.disturbanceSigma = [5; 5];              % Diffusion coefficients
p6.measurementCovariance = eye(2);         % Sensor covariance [cm^2]
p6.measurementVariance = 1;

sim6 = sim;
sim6.Ts = 10;
op6 = computeSteadyState([300; 300], p6.FMean, p6);
sim6.initialState = op6.x;
sim6.initialDisturbance = p6.FMean;
inputDeviation6 = zeros(2, sim6.N);        % Pumps stay at their nominal flows
measurementCovariance6 = eye(2);
P0Scale6 = 5;                              % Same documented demo initialization, P0 = 5 I

assert(abs(sim6.Ts / sim6.dt - round(sim6.Ts / sim6.dt)) < 1e-12, ...
    "Problem 6 sampling interval must be an integer number of SDE steps.");
sampleStride6 = round(sim6.Ts / sim6.dt);
sampleIndex6 = 1:sampleStride6:round((sim6.tf - sim6.t0) / sim6.dt) + 1;
assert(numel(sampleIndex6) == sim6.N, "SDE and KF time grids do not align.");
sampleTimes6 = sim6.t;

%% 6.2 P4 10% identified model and four-state realization
% The gains are the 10% identified values reported in Chapter 5. The shared
% time constants are the 10% Problem 4 values reported in Chapter 4.
K11 = 0.1336;
K12 = 0.0659;
K21 = 0.1181;
K22 = 0.2087;
T1 = 85.63;
T2 = 114.47;
T3 = 35.84;
T4 = 56.78;

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
[A4z, B4z] = c2dzoh(A4c, B4c, sim6.Ts);

% Chapter 6 Hankel/SVD realization. 89 Markov parameters are generated as reported.
% For a square N=44 Hankel, the unshifted/shifted blocks use H1:H87/H2:H88;
% the additional reported terms are retained for the full realization check.
nMarkovP4 = 89;
H4transfer = computeMarkovParameters(A4z, B4z, C4c, nMarkovP4);
hankelSizeP4 = 44;
HankelP4 = zeros(2 * hankelSizeP4, 2 * hankelSizeP4);
HankelShiftP4 = zeros(size(HankelP4));
for blockRow = 1:hankelSizeP4
    for blockColumn = 1:hankelSizeP4
        markovIndex = blockRow + blockColumn - 1;
        rows = (2 * blockRow - 1):(2 * blockRow);
        columns = (2 * blockColumn - 1):(2 * blockColumn);
        HankelP4(rows, columns) = H4transfer(:, :, markovIndex + 1);
        HankelShiftP4(rows, columns) = H4transfer(:, :, markovIndex + 2);
    end
end
[U4, S4, V4] = svd(HankelP4, "econ");
singularValuesP4 = diag(S4);
assert(size(HankelP4, 1) == 88 && numel(singularValuesP4) >= 5, ...
    "P4 Hankel dimensions are inconsistent with N=44.");
assert(singularValuesP4(4) > 0 && singularValuesP4(5) / singularValuesP4(4) < 1e-8, ...
    "P4 Markov Hankel matrix does not have the expected four-state order.");
SigmaRootP4 = diag(sqrt(singularValuesP4(1:4)));
SigmaInvRootP4 = diag(1 ./ sqrt(singularValuesP4(1:4)));
observabilityP4 = U4(:, 1:4) * SigmaRootP4;
controllabilityP4 = SigmaRootP4 * V4(:, 1:4)';
A4id = SigmaInvRootP4 * U4(:, 1:4)' * HankelShiftP4 * V4(:, 1:4) * SigmaInvRootP4;
B4id = controllabilityP4(:, 1:2);
C4id = observabilityP4(1:2, :);
H4realized = computeMarkovParameters(A4id, B4id, C4id, nMarkovP4);
markovDifferenceP4 = H4realized(:, :, 2:end) - H4transfer(:, :, 2:end);
markovTargetP4 = H4transfer(:, :, 2:end);
markovRelativeErrorP4 = norm(markovDifferenceP4(:)) / max(norm(markovTargetP4(:)), eps);
assert(size(A4id, 1) == 4 && markovRelativeErrorP4 < 1e-7, ...
    "P4 four-state realization failed its Markov-sequence check.");
fprintf("P4 Hankel order = 4; Markov relative error = %.2e.\n", markovRelativeErrorP4);

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

lin6 = computeJacobians(op6, p6);
aug6 = computeAugmentedModel(lin6, p6);
disc6 = discretizeModels(lin6, aug6, p6, sim6.Ts);
assert(max(abs(A6bar(:) - disc6.sde.A(:))) < 5.1e-4, ...
    "Copied P5 A matrix disagrees with computeAugmentedModel/discretizeModels.");
assert(max(abs(B6bar(:) - disc6.sde.B(:))) < 5.1e-4, ...
    "Copied P5 B matrix disagrees with computeAugmentedModel/discretizeModels.");
assert(max(abs(Q6bar(:) - disc6.sde.Q(:))) < 5.1e-2, ...
    "Copied P5 Q matrix disagrees with computeAugmentedModel/discretizeModels.");
assertCovarianceP6(Q6bar, "P5 Qbar");

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

%% 6.4 Problem 4 linear model, no step
% Simulate with the same Q and R used to design the filters.
rng(sim.randomSeed + 610);
modelP4Truth = modelP4;
modelP4Truth.P0 = zeros(4);
truthP4LinearNoStep = simulateLinearSystem(modelP4Truth, inputDeviation6, sim6.N);
resultP4LinearNoStep.dynamic = kalmanFilter( ...
    modelP4, truthP4LinearNoStep.y, inputDeviation6, "dynamic");
resultP4LinearNoStep.static = kalmanFilter( ...
    modelP4, truthP4LinearNoStep.y, inputDeviation6, "static");
p6Results.p4.linear.noStep.truth = truthP4LinearNoStep;
p6Results.p4.linear.noStep.filters = resultP4LinearNoStep;
p6Results.p4.linear.noStep.outputRMSE = [ ...
    sqrt(mean((modelP4.C * resultP4LinearNoStep.static.xf - ...
        modelP4.C * truthP4LinearNoStep.x).^2, "all")), ...
    sqrt(mean((modelP4.C * resultP4LinearNoStep.dynamic.xf - ...
        modelP4.C * truthP4LinearNoStep.x).^2, "all"))];
plotP6OutputFit(sampleTimes6, modelP4.C * truthP4LinearNoStep.x, ...
    truthP4LinearNoStep.y, modelP4.C, resultP4LinearNoStep.static, ...
    resultP4LinearNoStep.dynamic, "P4 linear no-step output estimates", ...
    fullfile(projectRoot, "figures", "problem6_p4_linear_no_step_outputs.png"));

%% 6.5 Problem 5 linear model, no step
rng(sim.randomSeed + 611);
modelP5Truth = modelP5;
modelP5Truth.P0 = zeros(6);
truthP5LinearNoStep = simulateLinearSystem(modelP5Truth, inputDeviation6, sim6.N);
resultP5LinearNoStep.dynamic = kalmanFilter( ...
    modelP5, truthP5LinearNoStep.y, inputDeviation6, "dynamic");
resultP5LinearNoStep.static = kalmanFilter( ...
    modelP5, truthP5LinearNoStep.y, inputDeviation6, "static");
p6Results.p5.linear.noStep.truth = truthP5LinearNoStep;
p6Results.p5.linear.noStep.filters = resultP5LinearNoStep;
p6Results.p5.linear.noStep.outputRMSE = [ ...
    sqrt(mean((modelP5.C * resultP5LinearNoStep.static.xf - ...
        modelP5.C * truthP5LinearNoStep.x).^2, "all")), ...
    sqrt(mean((modelP5.C * resultP5LinearNoStep.dynamic.xf - ...
        modelP5.C * truthP5LinearNoStep.x).^2, "all"))];
plotP6OutputFit(sampleTimes6, modelP5.C * truthP5LinearNoStep.x, ...
    truthP5LinearNoStep.y, modelP5.C, resultP5LinearNoStep.static, ...
    resultP5LinearNoStep.dynamic, "P5 linear no-step output estimates", ...
    fullfile(projectRoot, "figures", "problem6_p5_linear_no_step_outputs.png"));

%% 6.6 Nonlinear model, no step
rng(sim.randomSeed);
truthNoStepFine = simulateSDE_rerun(p6, sim6);
truthNoStep6 = truthNoStepFine.x(sampleIndex6, :);
truthNoStepD6 = truthNoStepFine.d(sampleIndex6, :);
measurementNoStep6 = truthNoStepFine.y(sampleIndex6, :).'- op6.y;
trueTankNoStep6 = truthNoStep6(:, 1:4).' - op6.x;
trueDisturbanceNoStep6 = truthNoStepD6.' - p6.FMean;
trueOutputNoStep6 = truthNoStep6(:, 1:2).' ./ (p6.rho .* p6.A(1:2)) - op6.y;
resultP4NoStep.dynamic = kalmanFilter(modelP4, measurementNoStep6, inputDeviation6, "dynamic");
resultP4NoStep.static = kalmanFilter(modelP4, measurementNoStep6, inputDeviation6, "static");
resultP5NoStep.dynamic = kalmanFilter(modelP5, measurementNoStep6, inputDeviation6, "dynamic");
resultP5NoStep.static = kalmanFilter(modelP5, measurementNoStep6, inputDeviation6, "static");
p6Results.nonlinear.noStep.truth = truthNoStepFine;
p6Results.nonlinear.noStep.filters.p4 = resultP4NoStep;
p6Results.nonlinear.noStep.filters.p5 = resultP5NoStep;

assert(isequal(measurementNoStep6, truthNoStepFine.y(sampleIndex6, :).'-op6.y) && ...
    isequal(inputDeviation6, zeros(size(inputDeviation6))), ...
    "The filters must use the same sampled measurements and known input deviations.");
assert(all(isfinite(resultP4NoStep.dynamic.xf(:))) && all(isfinite(resultP5NoStep.dynamic.xf(:))), ...
    "No-step filter returned nonfinite state estimates.");

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
assert(gainRelativeErrorP4 < 1e-5 && gainRelativeErrorP5 < 1e-5, ...
    "Dynamic no-step filter gains did not converge to their static gains.");

fprintf("\nProblem 6 no-step results (RMSE in native state/output units):\n");
fprintf("P4 abstract state: dynamic/static estimate difference RMSE %.4g; gain gap %.2e.\n", ...
    sqrt(mean((resultP4NoStep.dynamic.xf - resultP4NoStep.static.xf).^2, "all")), gainRelativeErrorP4);
fprintf("P5 physical state: dynamic/static estimate difference RMSE %.4g; gain gap %.2e.\n", ...
    sqrt(mean((resultP5NoStep.dynamic.xf - resultP5NoStep.static.xf).^2, "all")), gainRelativeErrorP5);
fprintf("P4 output vs physical truth: dynamic %.4g cm, static %.4g cm.\n", ...
    sqrt(mean((modelP4.C * resultP4NoStep.dynamic.xf - trueOutputNoStep6).^2, "all")), ...
    sqrt(mean((modelP4.C * resultP4NoStep.static.xf - trueOutputNoStep6).^2, "all")));
fprintf("P5 tank deviations: dynamic %.4g g, static %.4g g; D deviations: dynamic %.4g, static %.4g cm^3/s.\n", ...
    sqrt(mean((resultP5NoStep.dynamic.xf(1:4, :) - trueTankNoStep6).^2, "all")), ...
    sqrt(mean((resultP5NoStep.static.xf(1:4, :) - trueTankNoStep6).^2, "all")), ...
    sqrt(mean((resultP5NoStep.dynamic.xf(5:6, :) - trueDisturbanceNoStep6).^2, "all")), ...
    sqrt(mean((resultP5NoStep.static.xf(5:6, :) - trueDisturbanceNoStep6).^2, "all")));
fprintf("P5 output vs physical truth: dynamic %.4g cm, static %.4g cm; gain gap %.2e.\n", ...
    sqrt(mean((modelP5.C * resultP5NoStep.dynamic.xf - trueOutputNoStep6).^2, "all")), ...
    sqrt(mean((modelP5.C * resultP5NoStep.static.xf - trueOutputNoStep6).^2, "all")), gainRelativeErrorP5);

plotP6EstimatePair(sampleTimes6, resultP4NoStep.static, resultP4NoStep.dynamic, ...
    1:4, ["Abstract state 1", "Abstract state 2", "Abstract state 3", "Abstract state 4"], ...
    "P4 no-step: static vs dynamic estimates (not physical tank states)", ...
    fullfile(projectRoot, "figures", "problem6_p4_no_step_estimates.png"));
plotP6EstimatePair(sampleTimes6, resultP5NoStep.static, resultP5NoStep.dynamic, ...
    1:4, ["X_1", "X_2", "X_3", "X_4"], "P5 no-step: tank deviation estimates", ...
    fullfile(projectRoot, "figures", "problem6_p5_no_step_tank_estimates.png"));
plotP6TruthPair(sampleTimes6, trueTankNoStep6, resultP5NoStep.static.xf(1:4, :), ...
    resultP5NoStep.dynamic.xf(1:4, :), ["X_1", "X_2", "X_3", "X_4"], ...
    "P5 no-step: tank deviations vs physical truth", ...
    fullfile(projectRoot, "figures", "problem6_p5_no_step_tank_truth.png"));
plotP6EstimatePair(sampleTimes6, resultP5NoStep.static, resultP5NoStep.dynamic, ...
    5:6, ["D_3", "D_4"], "P5 no-step: disturbance deviation estimates", ...
    fullfile(projectRoot, "figures", "problem6_p5_no_step_disturbance_estimates.png"));
plotP6TruthPair(sampleTimes6, trueDisturbanceNoStep6, resultP5NoStep.static.xf(5:6, :), ...
    resultP5NoStep.dynamic.xf(5:6, :), ["D_3", "D_4"], ...
    "P5 no-step: disturbance deviations vs physical truth", ...
    fullfile(projectRoot, "figures", "problem6_p5_no_step_disturbance_truth.png"));

figure("Visible", "off");
semilogy(sampleTimes6, max(gainGapP4, eps), "LineWidth", 1.2);
hold on;
semilogy(sampleTimes6, max(gainGapP5, eps), "LineWidth", 1.2);
grid on;
xlabel("Time [s]");
ylabel("||K_{dynamic} - K_{static}||_F");
legend("P4", "P5", "Location", "best");
title("No-step dynamic Kalman gains converge to static gains");
saveas(gcf, fullfile(projectRoot, "figures", "problem6_no_step_gain_convergence.png"));
close(gcf);

plotP6OutputFit(sampleTimes6, trueOutputNoStep6, measurementNoStep6, modelP4.C, ...
    resultP4NoStep.static, resultP4NoStep.dynamic, "P4 no-step output fit", ...
    fullfile(projectRoot, "figures", "problem6_p4_no_step_outputs.png"));
plotP6OutputFit(sampleTimes6, trueOutputNoStep6, measurementNoStep6, modelP5.C, ...
    resultP5NoStep.static, resultP5NoStep.dynamic, "P5 no-step output fit", ...
    fullfile(projectRoot, "figures", "problem6_p5_no_step_outputs.png"));

%% 6.7 Offset-free filter models and step setup
stepTime6 = 12 * 60;
stepMean6 = [80; 100];
meanSchedule6.times = [sim6.t0; stepTime6];
meanSchedule6.values = [p6.FMean, stepMean6];
rng(sim.randomSeed);
truthStepFine = simulateSDE_rerun(p6, sim6, meanSchedule6);
truthStep6 = truthStepFine.x(sampleIndex6, :);
truthStepD6 = truthStepFine.d(sampleIndex6, :);
measurementStep6 = truthStepFine.y(sampleIndex6, :).'- op6.y;
trueTankStep6 = truthStep6(:, 1:4).' - op6.x;
trueDisturbanceStep6 = truthStepD6.' - p6.FMean;
trueOutputStep6 = truthStep6(:, 1:2).' ./ (p6.rho .* p6.A(1:2)) - op6.y;
trueMeanShift6 = zeros(2, sim6.N);
trueMeanShift6(:, sampleTimes6 >= stepTime6) = repmat(stepMean6 - p6.FMean, ...
    1, sum(sampleTimes6 >= stepTime6));
stepSampleIndex6 = find(sampleTimes6 == stepTime6, 1);
assert(~isempty(stepSampleIndex6) && ...
    isequal(truthNoStepFine.v(sampleIndex6, :), truthStepFine.v(sampleIndex6, :)), ...
    "No-step and stepped experiments must use matching seeded measurement-noise samples.");
assert(max(abs(truthNoStepFine.t(sampleIndex6) - sampleTimes6)) < 1e-10, ...
    "SDE samples must align with the 10-second KF grid.");
assert(norm(truthStepD6(stepSampleIndex6, :) - truthNoStepD6(stepSampleIndex6, :)) < 1e-12 && ...
    norm(truthStepD6(stepSampleIndex6 + 1, :) - truthNoStepD6(stepSampleIndex6 + 1, :)) > 0, ...
    "The disturbance mean change must start at exactly 720 seconds.");

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
SigmaNValues6 = [0.1, 1, 5];              % Per-sample standard deviations [cm^3/s]
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
    assert(norm(modelP5N{noiseIndex}.A(7:8, 7:8) - eye(2), "fro") < 1e-12, ...
        "P5 offset states must be discrete random walks.");
end

%% 6.8 Problem 4 and 5 linear models, step
stepTransient = sampleTimes6 >= stepTime6 & sampleTimes6 < stepTime6 + 180;
stepSettled = sampleTimes6 >= stepTime6 + 180;
postStepSummary6 = zeros(numel(SigmaNValues6), 12);
linearStepSummary6 = zeros(numel(SigmaNValues6), 5);
for noiseIndex = 1:numel(SigmaNValues6)
    sigmaN = SigmaNValues6(noiseIndex);
    modelP4N.Q = blkdiag(modelP4.Q, sigmaN^2 * eye(2));
    modelP4N.P0 = P0Scale6 * eye(6);
    assertKalmanModelP6(modelP4N);

    % The P4 step is an effective input bias, not a physical F3/F4 estimate.
    % A one-sample impulse sets the augmented bias state at the step.
    stepImpulse6 = zeros(2, sim6.N);
    stepImpulse6(:, stepSampleIndex6 - 1) = stepMean6 - p6.FMean;

    rng(sim.randomSeed + 620);
    modelP4NTruth = modelP4N;
    modelP4NTruth.P0 = zeros(6);
    modelP4NTruth.B = [modelP4N.B, [zeros(4, 2); eye(2)]];
    inputP4NTruth = [inputDeviation6; stepImpulse6];
    truthP4LinearStep = simulateLinearSystem( ...
        modelP4NTruth, inputP4NTruth, sim6.N);
    p4LinearStep.dynamic = kalmanFilter( ...
        modelP4N, truthP4LinearStep.y, inputDeviation6, "dynamic");
    p4LinearStep.static = kalmanFilter( ...
        modelP4N, truthP4LinearStep.y, inputDeviation6, "static");

    rng(sim.randomSeed + 621);
    modelP5NTruth = modelP5N{noiseIndex};
    modelP5NTruth.P0 = zeros(8);
    modelP5NTruth.B = [modelP5NTruth.B, [zeros(6, 2); eye(2)]];
    inputP5NTruth = [inputDeviation6; stepImpulse6];
    truthP5LinearStep = simulateLinearSystem( ...
        modelP5NTruth, inputP5NTruth, sim6.N);
    p5LinearStep.dynamic = kalmanFilter( ...
        modelP5N{noiseIndex}, truthP5LinearStep.y, inputDeviation6, "dynamic");
    p5LinearStep.static = kalmanFilter( ...
        modelP5N{noiseIndex}, truthP5LinearStep.y, inputDeviation6, "static");

    p6Results.p4.linear.step{noiseIndex}.truth = truthP4LinearStep;
    p6Results.p4.linear.step{noiseIndex}.filters = p4LinearStep;
    p6Results.p5.linear.step{noiseIndex}.truth = truthP5LinearStep;
    p6Results.p5.linear.step{noiseIndex}.filters = p5LinearStep;
    p6Results.p4.linear.step{noiseIndex}.sigmaN = sigmaN;
    p6Results.p5.linear.step{noiseIndex}.sigmaN = sigmaN;

    p4LinearOutput = modelP4N.C * truthP4LinearStep.x;
    p5LinearOutput = modelP5N{noiseIndex}.C * truthP5LinearStep.x;
    p5LinearStepRMSE = [ ...
        sqrt(mean((modelP5N{noiseIndex}.C * p5LinearStep.static.xp(:, stepTransient) - ...
            p5LinearOutput(:, stepTransient)).^2, "all")), ...
        sqrt(mean((modelP5N{noiseIndex}.C * p5LinearStep.dynamic.xp(:, stepTransient) - ...
            p5LinearOutput(:, stepTransient)).^2, "all"))];
    p4LinearStepRMSE = [ ...
        sqrt(mean((modelP4N.C * p4LinearStep.static.xp(:, stepTransient) - ...
            p4LinearOutput(:, stepTransient)).^2, "all")), ...
        sqrt(mean((modelP4N.C * p4LinearStep.dynamic.xp(:, stepTransient) - ...
            p4LinearOutput(:, stepTransient)).^2, "all"))];
    linearStepSummary6(noiseIndex, :) = [sigmaN, p4LinearStepRMSE, p5LinearStepRMSE];
    p6Results.p4.linear.step{noiseIndex}.transientOutputRMSE = p4LinearStepRMSE;
    p6Results.p5.linear.step{noiseIndex}.transientOutputRMSE = p5LinearStepRMSE;

    sigmaTag = strrep(num2str(sigmaN), ".", "p");
    plotP6OutputFit(sampleTimes6, p4LinearOutput, truthP4LinearStep.y, ...
        modelP4N.C, p4LinearStep.static, p4LinearStep.dynamic, ...
        "P4 linear step output estimates", fullfile(projectRoot, "figures", ...
        "problem6_p4_linear_step_outputs_sigma_" + sigmaTag + ".png"));
    plotP6OutputFit(sampleTimes6, p5LinearOutput, truthP5LinearStep.y, ...
        modelP5N{noiseIndex}.C, p5LinearStep.static, p5LinearStep.dynamic, ...
        "P5 linear step output estimates", fullfile(projectRoot, "figures", ...
        "problem6_p5_linear_step_outputs_sigma_" + sigmaTag + ".png"));

    %% Nonlinear model, step: evaluate both augmented filters on shared truth
    p4Step.dynamic = kalmanFilter(modelP4N, measurementStep6, inputDeviation6, "dynamic");
    p4Step.static = kalmanFilter(modelP4N, measurementStep6, inputDeviation6, "static");
    p5Step.dynamic = kalmanFilter(modelP5N{noiseIndex}, measurementStep6, inputDeviation6, "dynamic");
    p5Step.static = kalmanFilter(modelP5N{noiseIndex}, measurementStep6, inputDeviation6, "static");
    p6Results.nonlinear.step{noiseIndex}.truth = truthStepFine;
    p6Results.nonlinear.step{noiseIndex}.filters.p4 = p4Step;
    p6Results.nonlinear.step{noiseIndex}.filters.p5 = p5Step;
    p6Results.nonlinear.step{noiseIndex}.sigmaN = sigmaN;
    assert(all(isfinite(p4Step.dynamic.xf(:))) && all(isfinite(p5Step.dynamic.xf(:))), ...
        "Stepped filter returned nonfinite state estimates.");

    [~, staticGainP4N] = kalmanStatic(modelP4N);
    [~, staticGainP5N] = kalmanStatic(modelP5N{noiseIndex});
    % Gains depend only on the model covariances, so zero dummy data isolate convergence.
    gainProbeSamples6 = 3000;
    gainProbeInput6 = zeros(2, gainProbeSamples6);
    gainProbeOutput6 = zeros(2, gainProbeSamples6);
    gainProbeP4N = kalmanFilter(modelP4N, gainProbeOutput6, gainProbeInput6, "dynamic");
    gainProbeP5N = kalmanFilter(modelP5N{noiseIndex}, gainProbeOutput6, gainProbeInput6, "dynamic");
    probeGapP4N = norm(gainProbeP4N.K(:, :, end) - staticGainP4N, "fro") / ...
        max(norm(staticGainP4N, "fro"), eps);
    probeGapP5N = norm(gainProbeP5N.K(:, :, end) - staticGainP5N, "fro") / ...
        max(norm(staticGainP5N, "fro"), eps);
    assert(probeGapP4N < 1e-5 && probeGapP5N < 1e-5, ...
        "Dynamic augmented-model gains did not converge to their static gains.");
    p4StepGainGap = norm(p4Step.dynamic.K(:, :, end) - staticGainP4N, "fro") / ...
        max(norm(staticGainP4N, "fro"), eps);
    p5StepGainGap = norm(p5Step.dynamic.K(:, :, end) - staticGainP5N, "fro") / ...
        max(norm(staticGainP5N, "fro"), eps);

    p4PredictionDynamic = modelP4N.C * p4Step.dynamic.xp;
    p4PredictionStatic = modelP4N.C * p4Step.static.xp;
    p5PredictionDynamic = modelP5N{noiseIndex}.C * p5Step.dynamic.xp;
    p5PredictionStatic = modelP5N{noiseIndex}.C * p5Step.static.xp;
    postStepSummary6(noiseIndex, :) = [sigmaN, ...
        sqrt(mean((p5PredictionDynamic(:, stepTransient) - trueOutputStep6(:, stepTransient)).^2, "all")), ...
        max(abs(mean(p5PredictionDynamic(:, stepSettled) - trueOutputStep6(:, stepSettled), 2))), ...
        sqrt(mean((p5Step.dynamic.xf(7:8, stepSettled) - trueMeanShift6(:, stepSettled)).^2, "all")), ...
        sqrt(mean((p5PredictionStatic(:, stepTransient) - trueOutputStep6(:, stepTransient)).^2, "all")), ...
        max(abs(mean(p5PredictionStatic(:, stepSettled) - trueOutputStep6(:, stepSettled), 2))), ...
        sqrt(mean((p5Step.static.xf(7:8, stepSettled) - trueMeanShift6(:, stepSettled)).^2, "all")), ...
        sqrt(mean((p4PredictionDynamic(:, stepTransient) - trueOutputStep6(:, stepTransient)).^2, "all")), ...
        max(abs(mean(p4PredictionDynamic(:, stepSettled) - trueOutputStep6(:, stepSettled), 2))), ...
        sqrt(mean((p4PredictionStatic(:, stepTransient) - trueOutputStep6(:, stepTransient)).^2, "all")), ...
        max(abs(mean(p4PredictionStatic(:, stepSettled) - trueOutputStep6(:, stepSettled), 2))), ...
        max(p4StepGainGap, p5StepGainGap)];

    fprintf("\nStep offset sigma_N=%g cm^3/s:\n", sigmaN);
    fprintf("  P5 dynamic/static state difference RMSE: %.4g; dynamic gain gap: %.2e.\n", ...
        sqrt(mean((p5Step.dynamic.xf - p5Step.static.xf).^2, "all")), p5StepGainGap);
    fprintf("  P4 dynamic/static state difference RMSE: %.4g (bias states are effective pump biases); gain gap: %.2e.\n", ...
        sqrt(mean((p4Step.dynamic.xf - p4Step.static.xf).^2, "all")), p4StepGainGap);
    fprintf("  P5 dynamic transient output RMSE: %.4g cm; max abs settled output bias: %.4g cm; settled N RMSE: %.4g cm^3/s.\n", ...
        postStepSummary6(noiseIndex, 2:4));
    fprintf("  P5 static transient output RMSE: %.4g cm; max abs settled output bias: %.4g cm; settled N RMSE: %.4g cm^3/s.\n", ...
        postStepSummary6(noiseIndex, 5:7));
    fprintf("  P4 dynamic/static transient output RMSE: %.4g / %.4g cm; max abs settled output bias: %.4g / %.4g cm.\n", ...
        postStepSummary6(noiseIndex, 8), postStepSummary6(noiseIndex, 10), ...
        postStepSummary6(noiseIndex, 9), postStepSummary6(noiseIndex, 11));
    fprintf("  P5 tank RMSE dynamic/static: %.4g / %.4g g; D RMSE: %.4g / %.4g cm^3/s; N RMSE: %.4g / %.4g cm^3/s.\n", ...
        sqrt(mean((p5Step.dynamic.xf(1:4, :) - trueTankStep6).^2, "all")), ...
        sqrt(mean((p5Step.static.xf(1:4, :) - trueTankStep6).^2, "all")), ...
        sqrt(mean((p5Step.dynamic.xf(5:6, :) - trueDisturbanceStep6).^2, "all")), ...
        sqrt(mean((p5Step.static.xf(5:6, :) - trueDisturbanceStep6).^2, "all")), ...
        sqrt(mean((p5Step.dynamic.xf(7:8, :) - trueMeanShift6).^2, "all")), ...
        sqrt(mean((p5Step.static.xf(7:8, :) - trueMeanShift6).^2, "all")));

    plotP6EstimatePair(sampleTimes6, p4Step.static, p4Step.dynamic, 1:6, ...
        ["Abstract state 1", "Abstract state 2", "Abstract state 3", "Abstract state 4", ...
         "Effective input bias 1", "Effective input bias 2"], ...
        "P4 stepped: static vs dynamic (abstract states and pump biases)", ...
        fullfile(projectRoot, "figures", "problem6_p4_step_estimates_sigma_" + sigmaTag + ".png"));
    plotP6OutputFit(sampleTimes6, trueOutputStep6, measurementStep6, modelP4N.C, ...
        p4Step.static, p4Step.dynamic, "P4 stepped output fit", ...
        fullfile(projectRoot, "figures", "problem6_p4_step_outputs_sigma_" + sigmaTag + ".png"));

    plotP6EstimatePair(sampleTimes6, p5Step.static, p5Step.dynamic, 1:8, ...
        ["X_1", "X_2", "X_3", "X_4", "D_3", "D_4", "N_3", "N_4"], ...
        "P5 stepped: static vs dynamic estimates", ...
        fullfile(projectRoot, "figures", "problem6_p5_step_estimates_sigma_" + sigmaTag + ".png"));
    plotP6TruthPair(sampleTimes6, trueTankStep6, p5Step.static.xf(1:4, :), ...
        p5Step.dynamic.xf(1:4, :), ["X_1", "X_2", "X_3", "X_4"], ...
        "P5 stepped: tank deviations vs physical truth", ...
        fullfile(projectRoot, "figures", "problem6_p5_step_tank_truth_sigma_" + sigmaTag + ".png"));
    plotP6TruthPair(sampleTimes6, trueDisturbanceStep6, p5Step.static.xf(5:6, :), ...
        p5Step.dynamic.xf(5:6, :), ["D_3", "D_4"], ...
        "P5 stepped: disturbance deviations vs physical truth", ...
        fullfile(projectRoot, "figures", "problem6_p5_step_disturbance_truth_sigma_" + sigmaTag + ".png"));
    plotP6TruthPair(sampleTimes6, trueMeanShift6, p5Step.static.xf(7:8, :), ...
        p5Step.dynamic.xf(7:8, :), ["N_3", "N_4"], ...
        "P5 stepped: true and estimated mean shift", ...
        fullfile(projectRoot, "figures", "problem6_p5_step_mean_shift_sigma_" + sigmaTag + ".png"));
    plotP6OutputFit(sampleTimes6, trueOutputStep6, measurementStep6, modelP5N{noiseIndex}.C, ...
        p5Step.static, p5Step.dynamic, "P5 stepped output fit", ...
        fullfile(projectRoot, "figures", "problem6_p5_step_outputs_sigma_" + sigmaTag + ".png"));
end

fprintf("\nPost-step summary (first 180 s is transient; remaining interval is settled):\n");
disp(array2table(postStepSummary6, "VariableNames", ...
    ["sigmaN_cm3s", "P5_dynamic_transient_RMSE_cm", "P5_dynamic_maxAbsSettledBias_cm", ...
     "P5_dynamic_settled_N_RMSE_cm3s", "P5_static_transient_RMSE_cm", ...
     "P5_static_maxAbsSettledBias_cm", "P5_static_settled_N_RMSE_cm3s", ...
     "P4_dynamic_transient_RMSE_cm", "P4_dynamic_maxAbsSettledBias_cm", ...
     "P4_static_transient_RMSE_cm", "P4_static_maxAbsSettledBias_cm", ...
     "max_dynamic_static_gain_gap"]));
fprintf("\nLinear-model step output RMSE (columns: sigma_N, P4 static/dynamic, P5 static/dynamic):\n");
disp(linearStepSummary6);

%% 6.9 Compare and plot output estimation errors
nonlinearNoStepRMSE6 = [ ...
    sqrt(mean((modelP4.C * resultP4NoStep.static.xf - trueOutputNoStep6).^2, "all")), ...
    sqrt(mean((modelP4.C * resultP4NoStep.dynamic.xf - trueOutputNoStep6).^2, "all")); ...
    sqrt(mean((modelP5.C * resultP5NoStep.static.xf - trueOutputNoStep6).^2, "all")), ...
    sqrt(mean((modelP5.C * resultP5NoStep.dynamic.xf - trueOutputNoStep6).^2, "all"))];
noStepRMSEComparison6 = [p6Results.p4.linear.noStep.outputRMSE; ...
    p6Results.p5.linear.noStep.outputRMSE; nonlinearNoStepRMSE6];

middleNoiseIndex6 = find(SigmaNValues6 == 1, 1);
stepRMSEComparison6 = [linearStepSummary6(middleNoiseIndex6, 2:3); ...
    linearStepSummary6(middleNoiseIndex6, 4:5); ...
    postStepSummary6(middleNoiseIndex6, [10, 8]); ...
    postStepSummary6(middleNoiseIndex6, [5, 2])];
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
legend("Static", "Dynamic", "Location", "best");
title("Step output estimation, \sigma_N = 1");
exportgraphics(gcf, fullfile(projectRoot, "figures", ...
    "problem6_output_rmse_comparison.png"), "Resolution", 400);
close(gcf);