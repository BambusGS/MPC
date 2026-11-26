%% Problem 8: Unconstrained MPC
% Run problem4_main.m and then main_problem5.m first; they create G_final and disc in the workspace. This keeps model identification outside this script

root = fileparts(mfilename("fullpath"));
addpath(fullfile(root, "functions"));

% if ~exist("G_final", "var") || ~exist("disc", "var")
%     error("problem8:missingModels", "Run problem4_main.m and main_problem5.m before problem8.m.");
% end

% Convert the identified transfer functions to the same discrete state-space form used by the MPC functions. The chapter's prediction model has no D term
p = createParameters();
sim = createSimSettings(p);
models = struct("name", {}, "A", {}, "B", {}, "C", {}, "G", {});
stepPercent = [10, 25, 50];
if numel(G_final) ~= numel(stepPercent)
    error("problem8:identifiedModelCount", "Expected the three identified Problem 4 models.");
end

for i = 1:numel(G_final)
    [A, B, C, D] = ssdata(c2d(ss(G_final{i}), sim.Ts, "zoh"));
    if norm(D, "fro") > 1e-9
        error("problem8:directFeedthrough", "The MPC prediction model assumes zero direct feedthrough.");
    end
    models(end + 1) = struct("name", sprintf("Problem 4, %.0f%% step model", stepPercent(i)), "A", A, "B", B, "C", C, "G", zeros(size(A, 1), 1));
end

% % Include the discrete deterministic linearized model from Problem 5
% if abs(disc.Ts - sim.Ts) > 1e-9
%     error("problem8:samplingTime", "The Problem 5 model sampling time must match sim.Ts.");
% end
% models(end + 1) = struct("name", "Problem 5, deterministic linearized model", "A", disc.det.A, "B", disc.det.B, "C", disc.det.C, "G", disc.det.E);


%% Sim parameters setup
predictionTimeMinutes = 5;
requestedSteps = predictionTimeMinutes * 60 / sim.Ts;
predictionHorizon = round(requestedSteps);
actualPredictionTimeMinutes = predictionHorizon * sim.Ts / 60;

if predictionHorizon < 1
    error("problem8:predictionHorizon", "Prediction time must correspond to at least one sample.");
end
fprintf("Prediction horizon: %.2f min requested; %d steps = %.2f min actual (Ts = %.1f s).\n", predictionTimeMinutes, predictionHorizon, actualPredictionTimeMinutes, sim.Ts);


% Use deviation variables and a small output step as a simple closed-loop demonstration. Replace xhat and what with Kalman-filter estimates in use
simulationSteps = numel(sim.t) - 1;
stepTime = 15 * 60; % seconds
results = struct("name", {}, "time", {}, "output", {}, "input", {});

for i = 1:numel(models)
    model = models(i);
    n = size(model.A, 1);
    m = size(model.B, 2);
    pOut = size(model.C, 1);
    q = size(model.G, 2);

    Q = eye(pOut);
    S = 1e-2 * eye(m); % 100x smaller than the output weight (should not matter anyway here)
    ctrl = mpc_design(model.A, model.B, model.C, model.G, predictionHorizon, Q, S);

    % The example assumes the full state is available and the disturbance estimate is zero; the compute function also accepts filtered estimates
    xhat = zeros(n, 1);
    what = zeros(q, 1);
    uPrev = zeros(m, 1);
    
    % Stack references for the next N predicted output samples
    referenceTrajectory = zeros(pOut, numel(sim.t));
    referenceTrajectory(1, sim.t >= stepTime) = 10; % Tank 1: +10 cm
    
    y = zeros(pOut, simulationSteps + 1);
    u = zeros(m, simulationSteps);
    y(:, 1) = model.C * xhat;
    
    for k = 1:simulationSteps
        % Make step step reference trajectory for the next N predicted output samples
        refIndices = min(k + (1:predictionHorizon), numel(sim.t));
        R = reshape(referenceTrajectory(:, refIndices), [], 1);

        [Ustar, uApply, info] = mpc_compute(ctrl, xhat, what, R, uPrev);
        [Xpred, Ypred] = mpc_predict(model.A, model.B, model.C, model.G, xhat, what, Ustar);
        
        % Check both the affine control law and recursive prediction against their batch formulations to catch indexing or sign mistakes
        uGain = ctrl.Lx * xhat + ctrl.Lw * what + ctrl.LR * R + ctrl.Lu * uPrev;
        if norm(uApply - uGain, inf) > 1e-9
            error("problem8:gainMismatch", "The computed input does not match the affine control law.");
        end
        if norm(Ypred(:) - (info.b + ctrl.Gamma * Ustar), inf) > 1e-9
            error("problem8:predictionMismatch", "Recursive and batch output predictions do not match.");
        end

        xhat = Xpred(:, 1);
        y(:, k + 1) = Ypred(:, 1);
        u(:, k) = uApply;
        uPrev = uApply;
    end

    nominalFlow = [sim.F1; sim.F2];
    appliedFlow = u + repmat(nominalFlow, 1, size(u, 2));

    finalError = referenceTrajectory(:, end) - y(:, end);
    fprintf("%s: final tracking error = %.4g cm\n", model.name, norm(finalError));


    % Plot the results for this model [TBD, this is quick and dirty]
    % Compare predicted outputs with their references and show held inputs.
    figure("Name", model.name);
    tiledlayout(pOut + m, 1);
    sgtitle(model.name);

    for j = 1:pOut
        nexttile;
        plot(sim.t / 60, y(j,:).', "LineWidth", 1.5);
        hold on;
        plot(sim.t / 60, referenceTrajectory(j,:).', "--", "LineWidth", 1.2);
        ylabel(sprintf("Tank %d level deviation [cm]", j));
        legend("Output", "Reference", "Location", "best");
        grid on;
    end

    % Extend the last input sample to the final time for the stairs plot.
    inputPlot = [appliedFlow, appliedFlow(:, end)];
    for j = 1:m
        nexttile;
        stairs(sim.t / 60, inputPlot(j,:).', "LineWidth", 1.5);
        ylabel(sprintf("Pump flow F_%d [cm^3/s]", j));
        grid on;
    end
    xlabel("Time [min]");

end