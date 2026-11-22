rootDir = fileparts(mfilename('fullpath'));
addpath(fullfile(rootDir, 'functions'), '-begin');   % our copies take precedence
figDir = fullfile(rootDir, 'figures');
if ~exist(figDir, 'dir'), mkdir(figDir); end

p = createParameters();
sim = createSimSettings(p);
Ts = sim.Ts;

%% 6.1 Discrete-time state space models

% Problem 5 model: operating point, linearization and SDE model (Chapter 5)
us = [sim.F1; sim.F2];
op = computeSteadyState(us, p.FMean, p);
lin = computeJacobians(op, p);
aug = computeAugmentedModel(lin, p);

% Problem 4 model: identified 10 % fit with shared T1 and T2
if exist('estimated_params', 'var')
    ident = sharedTimeConstants(estimated_params, 1);
    disp('Problem 4 model: 10 % fit from problem4_main.m (workspace).');
else
    ident.K = [0.1336, 0.0659; 0.1181, 0.2087];             % Table 5.3, 10 %
    ident.T = [mean([85.6, 85.5]); mean([115.4, 114.5]); 36.0; 56.1];
    disp('Problem 4 model: 10 % values of Table 5.3 (run problem4_main.m first to use the fit).');
end
fprintf('  K = [%.4f %.4f; %.4f %.4f], T1..T4 = [%.1f %.1f %.1f %.1f] s\n', ...
    ident.K(1, 1), ident.K(1, 2), ident.K(2, 1), ident.K(2, 2), ident.T);

% Tuning (design parameters):
sigmaW4 = (p.disturbanceSigma ./ sqrt(2 * p.disturbanceRate)) ./ [1 - p.gamma(2); 1 - p.gamma(1)];
sigmaN5 = [0.3; 0.3];
sigmaN4 = [3; 3];

models = {
    buildKalmanModelP5(aug, lin, p, Ts, false)
    buildKalmanModelP5(aug, lin, p, Ts, true, sigmaN5)
    buildKalmanModelP4(ident, p, Ts, false, sigmaW4)
    buildKalmanModelP4(ident, p, Ts, true, sigmaW4, sigmaN4)
};

disp('Kalman filter models:');
for i = 1:numel(models)
    model = models{i};
    [~, ~, ~, ~, info] = kalmanStatic(model);
    fprintf('  %-28s n = %d | detectable = %d | max|eig(A - Kpx*C)| = %.4f | %s', ...
        model.name, size(model.A, 1), isDetectable(model.A, model.C), ...
        max(abs(info.eigErrorDynamics)), info.method);
    if strcmp(info.method, 'fixed-point iteration')
        fprintf(' (%d iterations)\n', info.iterations);
    else
        fprintf(' (rel. diff. to the Riccati iteration %.1e)\n', info.relDiffIteration);
    end
end

%% Scenario
tEnd = 60 * 60;                     % 60 min
t = 0:Ts:tEnd;
N = numel(t);

u = repmat(us, 1, N);               % pump steps to excite the system
u(1, t >= 10 * 60) = us(1) + 30;    % +10 % on pump 1 at 10 min
u(2, t >= 40 * 60) = us(2) - 30;    % -10 % on pump 2 at 40 min
U = u - us;                         % deviation inputs for the filters

FMeanNoStep = repmat(p.FMean, 1, N);
FMeanStep = FMeanNoStep;
FMeanStep(1, t >= 25 * 60) = p.FMean(1) + 30;   % step in the mean of F3 (as in Problem 3)

h0dev = [3; 3; 3; 3];               % plant starts 3 cm off the steady state,
                                    % the filters start at x_hat(0|-1) = 0
rng(sim.randomSeed);
noise = createKalmanNoise(p, Ts, N, 0.1);   % same noise for all runs

%% 6.2-6.4 Simulations
plantTypes = {'linear', 'nonlinear'};
scenarios = {'nostep', 'step'};
results = struct('scenario', {}, 'plant', {}, 'modelName', {}, 'model', {}, ...
    'dynamic', {}, 'static', {}, 'Xtrue', {});

for is = 1:numel(scenarios)
    if strcmp(scenarios{is}, 'step')
        FMean = FMeanStep;
        scenarioText = 'step in F_3 mean at 25 min';
    else
        FMean = FMeanNoStep;
        scenarioText = 'no steps';
    end
    for ip = 1:numel(plantTypes)
        plant = simulateKalmanPlant(p, op, lin, plantTypes{ip}, t, u, FMean, noise, h0dev);
        Y = plant.y - op.y;         % deviation measurements

        for im = 1:numel(models)
            model = models{im};
            if strcmp(scenarios{is}, 'nostep') && model.offsetFree
                continue            % 6.2: plain models only
            end
            resDyn = kalmanFilter(model, Y, U, 'dynamic');
            resStat = kalmanFilter(model, Y, U, 'static');
            Xtrue = truthInModelStates(p, op, plant, model);

            newResult.scenario = scenarios{is};
            newResult.plant = plantTypes{ip};
            newResult.modelName = model.name;
            newResult.model = model;
            newResult.dynamic = resDyn;
            newResult.static = resStat;
            newResult.Xtrue = Xtrue;
            results(end + 1) = newResult; %#ok<SAGROW>

            estimates = struct('x', {resDyn.xf, resStat.xf}, 'sd', {resDyn.sdf, []}, ...
                'label', {'Dynamic KF', 'Static KF'});
            plotKalmanStates(sprintf('%s - %s plant, %s', model.name, plantTypes{ip}, scenarioText), ...
                t, Xtrue, estimates, model, fullfile(figDir, sprintf('problem6_states_P%d_n%d_%s_%s.png', ...
                model.problem, size(model.A, 1), scenarios{is}, plantTypes{ip})));
        end
    end
end

%% Gain convergence (static vs dynamic), no step, linear plant
for im = [1, 3]
    model = models{im};
    r = findResult(results, 'nostep', 'linear', model.name);
    [~, Kstatic] = kalmanStatic(model);
    plotKalmanGains(sprintf('Gain convergence - %s', model.name), t, r.dynamic, Kstatic, model, ...
        fullfile(figDir, sprintf('problem6_gains_P%d_n%d.png', model.problem, size(model.A, 1))));
end

%% Innovations with the step: plain vs offset-free (dynamic filters)
pairs = {[1, 2], [3, 4]};
for ip = 1:numel(plantTypes)
    for ipair = 1:numel(pairs)
        res = cell(1, 2);
        labels = cell(1, 2);
        for j = 1:2
            model = models{pairs{ipair}(j)};
            r = findResult(results, 'step', plantTypes{ip}, model.name);
            res{j} = r.dynamic;
            labels{j} = model.name;
        end
        plotInnovations(sprintf('Innovations, %s plant, step in F_3 mean at 25 min', plantTypes{ip}), ...
            t, res, labels, fullfile(figDir, sprintf('problem6_innovations_P%d_%s.png', ...
            models{pairs{ipair}(1)}.problem, plantTypes{ip})));
    end
end

%% Results table
idxTransient = find(t <= 5 * 60);
idxAfter = find(t > 5 * 60);
idxLast = find(t > 40 * 60);
filterTypes = {'dynamic', 'static'};

fid = fopen(fullfile(figDir, 'problem6_results.csv'), 'w');
fprintf(fid, 'scenario,plant,model,filter,RMSE_transient_cm,RMSE_h1_cm,RMSE_h2_cm,offset_h1_cm,mean_e1_cm,mean_e2_cm,NIS\n');
fprintf('\n%-7s %-10s %-28s %-8s | %7s | %7s %7s | %8s | %7s %7s | %5s\n', 'Case', 'Plant', ...
    'Model', 'Filter', 'RMSE_tr', 'RMSE_h1', 'RMSE_h2', 'off_h1', 'mean_e1', 'mean_e2', 'NIS');
for i = 1:numel(results)
    r = results(i);
    for jf = 1:numel(filterTypes)
        res = r.(filterTypes{jf});
        err = r.model.scale(1:2) .* res.xf(1:2, :) - r.Xtrue(1:2, :);
        rmseTr = sqrt(mean(mean(err(:, idxTransient) .^ 2)));
        rmse = sqrt(mean(err(:, idxAfter) .^ 2, 2));
        offset = mean(err(1, idxLast));
        stats = innovationStats(res, idxAfter);
        fprintf('%-7s %-10s %-28s %-8s | %7.3f | %7.3f %7.3f | %8.3f | %7.3f %7.3f | %5.2f\n', ...
            r.scenario, r.plant, r.modelName, filterTypes{jf}, rmseTr, rmse, offset, stats.mean, stats.NIS);
        fprintf(fid, '%s,%s,%s,%s,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.3f\n', ...
            r.scenario, r.plant, r.modelName, filterTypes{jf}, rmseTr, rmse, offset, stats.mean, stats.NIS);
    end
end
fclose(fid);
fprintf('Figures and problem6_results.csv saved in %s\n', figDir);
