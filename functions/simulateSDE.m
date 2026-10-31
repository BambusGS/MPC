function result = simulateSDE(p, sim)
    % Simulate the augmented mean-reverting SDE

    % Simulate the augmented tank and disturbance states
    time = (sim.t0:sim.dt:sim.tf).';
    nSamples = numel(time);
    state = zeros(nSamples, 6);
    disturbance = zeros(nSamples, 2);
    state(1, 1:4) = sim.initialState.';
    disturbance(1, :) = sim.initialDisturbance.';
    pumpFlow = sim.input(:, 1);

    % Pre-compute the disturbance path with Euler, indeed Ts independent
    for k = 1:(nSamples - 1)
        currentDisturbance = disturbance(k, :).';
        dW = sqrt(sim.dt) .* randn(2, 1);
        disturbanceDerivative = p.disturbanceRate .* (p.FMean - currentDisturbance); % dF = a(F_bar - F)dt + sigma dW
        disturbance(k + 1, :) = (currentDisturbance + sim.dt .* disturbanceDerivative + p.disturbanceSigma .* dW).';
    end

    state(:, 5:6) = disturbance;

    % Integrate tank states with ode45 and calculate outputs
    measurementNoise = zeros(nSamples, 2);
    measuredOutput = zeros(nSamples, 2);
    systemOutputValues = zeros(nSamples, 2);
    for k = 1:nSamples  % each k is size Ts
        massState = state(k, 1:4).';
        disturbanceState = disturbance(k, :).';

        if k < nSamples
            [~, intervalState] = ode45(@(t, x) quadrupleTankProcess(t, x, pumpFlow, disturbanceState, p), [time(k), time(k + 1)], massState); % dx_m = f_det(x_m, u, x_d, p)dt
            state(k + 1, 1:4) = intervalState(end, :);
        end

        measurementNoise(k, :) = sqrt(diag(p.measurementCovariance)).' .* randn(1, 2); % v(t) ~ N(0, R)
        noiselessMeasurement = systemSense(massState, p); % g(x_m(t), p)
        measuredOutput(k, :) = (noiselessMeasurement + measurementNoise(k, :).').'; % g(x_m(t), p) + v(t)
        systemOutputValues(k, :) = systemOutput(massState, p).'; % z(t) = h(x_m(t), p)
    end

    % Store the simulation results
    result.t = time;
    result.x = state;
    result.d = state(:, 5:6);
    result.v = measurementNoise;
    result.y = measuredOutput;
    result.z = systemOutputValues;
end