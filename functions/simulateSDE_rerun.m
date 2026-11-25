function result = simulateSDE_rerun(p, sim, meanSchedule)
    % Rerun the nonlinear SDE with an optional piecewise-constant disturbance mean.

    if nargin < 3
        meanSchedule.times = sim.t0;
        meanSchedule.values = p.FMean;
    end

    time = (sim.t0:sim.dt:sim.tf).';
    nSamples = numel(time);
    state = zeros(nSamples, 6);
    disturbance = zeros(nSamples, 2);
    state(1, 1:4) = sim.initialState.';
    disturbance(1, :) = sim.initialDisturbance.';
    pumpFlow = sim.input(:, 1);

    meanIndex = 1;
    for k = 1:(nSamples - 1)
        currentDisturbance = disturbance(k, :).';
        dW = sqrt(sim.dt) .* randn(2, 1);
        while meanIndex < numel(meanSchedule.times) && time(k) >= meanSchedule.times(meanIndex + 1)
            meanIndex = meanIndex + 1;
        end
        currentMean = meanSchedule.values(:, meanIndex);
        disturbanceDerivative = p.disturbanceRate .* (currentMean - currentDisturbance);
        disturbance(k + 1, :) = (currentDisturbance + sim.dt .* disturbanceDerivative + ...
            p.disturbanceSigma .* dW).';
    end
    state(:, 5:6) = disturbance;

    measurementNoise = zeros(nSamples, 2);
    measuredOutput = zeros(nSamples, 2);
    systemOutputValues = zeros(nSamples, 2);
    for k = 1:nSamples
        massState = state(k, 1:4).';
        disturbanceState = disturbance(k, :).';
        if k < nSamples
            [~, intervalState] = ode45(@(t, x) quadrupleTankProcess(t, x, pumpFlow, disturbanceState, p), ...
                [time(k), time(k + 1)], massState);
            state(k + 1, 1:4) = intervalState(end, :);
        end

        measurementNoise(k, :) = sqrt(diag(p.measurementCovariance)).' .* randn(1, 2);
        noiselessMeasurement = systemSense(massState, p);
        measuredOutput(k, :) = (noiselessMeasurement + measurementNoise(k, :).').';
        systemOutputValues(k, :) = systemOutput(massState, p).';
    end

    result.t = time;
    result.x = state;
    result.d = state(:, 5:6);
    result.v = measurementNoise;
    result.y = measuredOutput;
    result.z = systemOutputValues;
end
