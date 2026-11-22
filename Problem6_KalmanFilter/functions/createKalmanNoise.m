function noise = createKalmanNoise(p, Ts, N, dt)
    % One realization of the Brownian increments and the measurement noise.
    % It is shared by the linear and the nonlinear plant so that both are
    % compared on the same disturbance trajectory. Call rng() before.
    %
    % Inputs:
    %   p: Parameter structure created by createParameters()
    %   Ts: Sampling time [s];  N: Number of samples;  dt: Fine time step [s]

    nSub = round(Ts / dt);
    noise.dt = dt;
    noise.dW = sqrt(dt) * randn(2, (N - 1) * nSub);
    noise.v = chol(p.measurementCovariance, 'lower') * randn(2, N);
end
