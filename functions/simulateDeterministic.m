function result = simulateDeterministic(p, sim)
    % Simulate the deterministic model with constant F3 and F4

    % Simulate the tank states
    time = sim.t;
    state = zeros(numel(time), 4); % numel gets the number of elements in the time vector
    state(1, :) = sim.initialState.';
    disturbance = sim.initialDisturbance;

    % Calculate measurements and outputs
    measuredOutput = zeros(numel(time), 2);
    systemOutputValues = zeros(numel(time), 2);
    for k = 1:numel(time) % each k is size Ts
        if k < numel(time)
            [~, intervalState] = ode45(@(t, x) quadrupleTankProcess(t, x, sim.input(:, k), disturbance, p), [time(k), time(k + 1)], state(k, :).'); % x_dot = f(x, u, d, p)
            state(k + 1, :) = intervalState(end, :);
        end

        measuredOutput(k, :) = systemSense(state(k, :).', p).'; % y(t) = g(x(t), p)
        systemOutputValues(k, :) = systemOutput(state(k, :).', p).'; % z(t) = h(x(t), p)
    end

    % Store the simulation results
    result.t = time;
    result.x = state;
    result.d = repmat(disturbance.', numel(time), 1);
    result.y = measuredOutput;
    result.z = systemOutputValues;
end