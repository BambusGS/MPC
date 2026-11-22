rootDir = fileparts(mfilename('fullpath'));
addpath(fullfile(rootDir, 'functions'), '-begin');
rng(1);

% Random model with correlated process and measurement noise
A = [0.90 0.10 0.00; 0.00 0.80 0.20; 0.10 0.00 0.70];
B = [1.0; 0.0; 0.5];
G = [1.0 0.0; 0.0 1.0; 0.5 0.5];
C = [1.0 0.0 0.0; 0.0 0.0 1.0];
L = [1.0 0 0 0; 0.3 0.8 0 0; 0.4 -0.2 0.7 0; 0.1 0.3 0.2 0.6];
Sigma = L * L';                      % joint covariance of [w; v]
model.A = A;  model.B = B;  model.G = G;  model.C = C;
model.Q = Sigma(1:2, 1:2);  model.S = Sigma(1:2, 3:4);  model.R = Sigma(3:4, 3:4);
model.P0 = 10 * eye(3);
model.name = 'random test model';

nPass = 0;  nTest = 0;

% 1. Riccati residual
[P, Kfx, Re, Pf, info] = kalmanStatic(model);
M = A * P * C' + G * model.S;
residual = norm(A * P * A' + G * model.Q * G' - M / (C * P * C' + model.R) * M' - P, 'fro') / norm(P, 'fro');
[nPass, nTest] = testCheck(residual < 1e-9, sprintf('1. DARE residual = %.1e', residual), nPass, nTest);

% 2. Solver vs iteration
[nPass, nTest] = testCheck(info.relDiffIteration < 1e-8, ...
    sprintf('2. %s vs fixed-point iteration: rel. diff. = %.1e', info.method, info.relDiffIteration), nPass, nTest);

% 3. Dynamic gain -> static gain
N = 3000;
x = zeros(3, N);  y = zeros(2, N);  u = sin((1:N) / 50);
xk = chol(model.P0, 'lower') * randn(3, 1);
for k = 1:N
    wv = L * randn(4, 1);
    x(:, k) = xk;
    y(:, k) = C * xk + wv(3:4);
    xk = A * xk + B * u(k) + G * wv(1:2);
end
resDyn = kalmanFilter(model, y, u, 'dynamic');
resStat = kalmanFilter(model, y, u, 'static');
gainDiff = norm(resDyn.K(:, :, end) - Kfx, 'fro') / norm(Kfx, 'fro');
[nPass, nTest] = testCheck(gainDiff < 1e-8, sprintf('3. ||Kfx,N - Kfx|| / ||Kfx|| = %.1e', gainDiff), nPass, nTest);
estDiff = max(max(abs(resDyn.xf(:, end-100:end) - resStat.xf(:, end-100:end))));
[nPass, nTest] = testCheck(estDiff < 1e-6, sprintf('   static = dynamic estimates after convergence: max diff = %.1e', estDiff), nPass, nTest);

% 4. S = 0: compare with the classical textbook filter
model0 = model;  model0.S = zeros(2, 2);
res0 = kalmanFilter(model0, y, u, 'dynamic');
xp = zeros(3, 1);  Pk = model0.P0;  maxDiff = 0;
for k = 1:N
    Kk = Pk * C' / (C * Pk * C' + model0.R);
    xf = xp + Kk * (y(:, k) - C * xp);
    Pfk = (eye(3) - Kk * C) * Pk;
    maxDiff = max(maxDiff, max(abs(xf - res0.xf(:, k))));
    xp = A * xf + B * u(k);
    Pk = A * Pfk * A' + G * model0.Q * G';
end
[nPass, nTest] = testCheck(maxDiff < 1e-8, sprintf('4. S = 0 vs textbook Kalman filter: max diff = %.1e', maxDiff), nPass, nTest);

% 5. Consistency (after the transient)
idx = 200:N;
stats = innovationStats(resDyn, idx);
err = x(:, idx) - resDyn.xf(:, idx);
covRatio = trace(cov(err')) / trace(Pf);
[nPass, nTest] = testCheck(abs(stats.NIS - 2) < 0.15, sprintf('5. NIS = %.3f (expected 2)', stats.NIS), nPass, nTest);
[nPass, nTest] = testCheck(abs(covRatio - 1) < 0.1, sprintf('   trace(cov(x - x_hat)) / trace(Pf) = %.3f (expected 1)', covRatio), nPass, nTest);
[nPass, nTest] = testCheck(all(abs(stats.rho1) < 0.06), sprintf('   lag-1 autocorrelation of e = [%.3f %.3f] (expected 0)', stats.rho1), nPass, nTest);

% 6. Four-tank models
p = createParameters();
sim = createSimSettings(p);
op = computeSteadyState([sim.F1; sim.F2], p.FMean, p);
lin = computeJacobians(op, p);
aug = computeAugmentedModel(lin, p);
ident.K = [0.1336, 0.0659; 0.1181, 0.2087];
ident.T = [85.55; 114.95; 36.0; 56.1];
tankModels = {buildKalmanModelP5(aug, lin, p, sim.Ts, false), buildKalmanModelP5(aug, lin, p, sim.Ts, true), ...
              buildKalmanModelP4(ident, p, sim.Ts, false), buildKalmanModelP4(ident, p, sim.Ts, true)};
for i = 1:numel(tankModels)
    m = tankModels{i};
    [~, ~, ~, ~, infoTank] = kalmanStatic(m);
    ok = isDetectable(m.A, m.C) && max(abs(infoTank.eigErrorDynamics)) < 1;
    [nPass, nTest] = testCheck(ok, sprintf('6. %s: detectable, max|eig(A - Kpx*C)| = %.4f', ...
        m.name, max(abs(infoTank.eigErrorDynamics))), nPass, nTest);
end

fprintf('\n%d of %d checks passed.\n', nPass, nTest);
if nPass < nTest
    error('testKalmanFilter:failed', 'Some checks failed.');
end
