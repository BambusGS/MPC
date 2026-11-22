function plant = simulateKalmanPlant(p, op, lin, plantType, t, u, FMean, noise, h0dev)
    % True plant for the Kalman filter tests, with SDE disturbances
    %
    % plantType = "linear":    linearized model (5.22), dX = A*X + B*U + E*D
    % plantType = "nonlinear": quadrupleTankProcess (2.10)
    % The disturbances follow the mean-reverting SDE of Section 2.3,
    %   dF = a_F*(FMean(t) - F)*dt + sigma_F*dw,
    % where FMean may contain steps (Problem 6.3). Euler-Maruyama on the fine
    % grid noise.dt; u and FMean are held constant over each sample (ZOH).
    %
    % Inputs:
    %   p: Parameter structure created by createParameters()
    %   op: Steady state created by computeSteadyState()
    %   lin: Linearized model created by computeJacobians()
    %   plantType: "linear" or "nonlinear"
    %   t: Sample times [s], 1 x N
    %   u: Pump flows [cm^3/s], 2 x N
    %   FMean: Mean disturbance flows [cm^3/s], 2 x N
    %   noise: Structure created by createKalmanNoise() (same for both plants)
    %   h0dev: Initial level deviations from the steady state [cm], 4 x 1
    %
    % Output:
    %   plant: Structure with the fields (at the sample times)
    %       h: Levels [cm], 4 x N;  F: Disturbance flows [cm^3/s], 2 x N
    %       q: Tank outflows [cm^3/s], 4 x N;  y: Measured levels [cm], 2 x N
    %       FMean, u, t, type

    if nargin < 9, h0dev = zeros(4, 1); end
    N = numel(t);
    Ts = t(2) - t(1);
    dt = noise.dt;
    nSub = round(Ts / dt);
    isLinear = strcmpi(plantType, 'linear');

    m = op.x + p.rho * p.A .* h0dev; % Tank masses [g]
    F = op.d;                        % Disturbance flows [cm^3/s]

    plant.h = zeros(4, N);
    plant.F = zeros(2, N);
    plant.q = zeros(4, N);
    plant.y = zeros(2, N);
    plant.FMean = FMean;
    plant.u = u;
    plant.t = t;
    plant.type = plantType;

    for k = 1:N
        h = max(m, 0) ./ (p.rho * p.A);
        plant.h(:, k) = h;
        plant.F(:, k) = F;
        if isLinear
            plant.q(:, k) = op.q + (m - op.x) ./ (p.rho * lin.T); % linearized outflows
        else
            plant.q(:, k) = p.a .* sqrt(2 * p.g * h);
        end
        plant.y(:, k) = h(1:2) + noise.v(:, k);
        if k == N
            break
        end

        dW = noise.dW(:, (k - 1) * nSub + (1:nSub));
        for j = 1:nSub
            if isLinear
                dm = lin.A * (m - op.x) + lin.B * (u(:, k) - op.u) + lin.E * (F - op.d);
            else
                dm = quadrupleTankProcess(0, m, u(:, k), F, p);
            end
            dF = p.disturbanceRate .* (FMean(:, k) - F) * dt + p.disturbanceSigma .* dW(:, j);
            m = m + dm * dt;
            F = F + dF;
            if ~isLinear
                m = max(m, 0);
                F = max(F, 0); % flows cannot be negative
            end
        end
    end
end
