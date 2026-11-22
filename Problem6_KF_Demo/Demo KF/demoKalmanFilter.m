%% Static and dynamic Kalman filter on a linear stochastic system
model.A = [0.90 0.10 0.00; 0.00 0.80 0.20; 0.10 0.00 0.70];
model.B = [1.0; 0.0; 0.5];
model.G = [1.0 0.0; 0.0 1.0; 0.5 0.5];
model.C = [1.0 0.0 0.0; 0.0 0.0 1.0];
model.Q = 0.5 * eye(2);
model.R = 0.2 * eye(2);
model.S = zeros(2, 2);
model.P0 = 5 * eye(3);

N = 300;
k = 0:N-1;
u = sin(k / 20);
rng(1);
sim = simulateLinearSystem(model, u, N);

resDyn = kalmanFilter(model, sim.y, u, 'dynamic');
resStat = kalmanFilter(model, sim.y, u, 'static');

plotStateEstimates(k, sim.x, {resDyn, resStat}, {'Dynamic KF', 'Static KF'}, ...
    {'x_1', 'x_2', 'x_3'}, 'Static and dynamic Kalman filter', 'demo_kalman_states.png');

[P, Kfx, Re, Pf, info] = kalmanStatic(model);
fprintf('Riccati equation solved with: %s\n', info.method);
fprintf('Static gain Kfx:\n'); disp(Kfx);
fprintf('||Kfx,N - Kfx|| / ||Kfx|| = %.1e\n', norm(resDyn.K(:, :, end) - Kfx, 'fro') / norm(Kfx, 'fro'));
stats = innovationStats(resDyn, 50:N);
fprintf('NIS = %.2f (expected %d)\n', stats.NIS, size(model.C, 1));
