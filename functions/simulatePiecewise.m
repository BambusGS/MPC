function result = simulatePiecewise(p, sim)
    % Simulate random F3/F4 values held constant between samples

    % Simulate the tank states
    time = sim.t;
    nSamples = numel(time);
    state = zeros(nSamples, 4);
    measurementNoise = zeros(nSamples, 2);
    disturbance = max(0, sim.initialDisturbance.' + p.piecewiseSigma.' .* randn(nSamples, 2));
    disturbance(1, :) = sim.initialDisturbance.';
    state(1, :) = sim.initialState.';

    measuredOutput = zeros(nSamples, 2);
    systemOutputValues = zeros(nSamples, 2);

    for k = 1:nSamples  % each k is size Ts
        if k < nSamples
            [~, intervalState] = ode45(@(t, x) quadrupleTankProcess(t, x, sim.input(:, k), disturbance(k, :).', p), [time(k), time(k + 1)], state(k, :).'); % x_dot = f(x, u, d, p)
            state(k + 1, :) = intervalState(end, :);
        end

        measurementNoise(k, :) = sqrt(diag(p.measurementCovariance)).' .* randn(1, 2); % v(t) ~ N(0, R)
        noiselessMeasurement = systemSense(state(k, :).', p); % g(x(t), p)
        measuredOutput(k, :) = (noiselessMeasurement + measurementNoise(k, :).').'; % g(x(t), p) + v(t)
        systemOutputValues(k, :) = systemOutput(state(k, :).', p).'; % z(t) = h(x(t), p)
    end

    % Store the simulation results
    result.t = time;
    result.x = state;
    result.d = disturbance;
    samplesPerInterval = round(sim.Ts / sim.dt);
    result.tPlot = (sim.t0:sim.dt:sim.tf).';
    result.dPlot = zeros(numel(result.tPlot), 2);
    for k = 1:(nSamples - 1)
        firstIndex = (k - 1) * samplesPerInterval + 1;
        lastIndex = k * samplesPerInterval;
        for j = firstIndex:lastIndex
            result.dPlot(j, :) = disturbance(k, :);
        end
    end
    result.dPlot(end, :) = disturbance(end, :);
    result.v = measurementNoise;
    result.y = measuredOutput;
    result.z = systemOutputValues;
end