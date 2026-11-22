function model = buildKalmanModelP5(aug, lin, p, Ts, offsetFree, sigmaN)
    % Discrete-time Kalman filter model from the Problem 5 SDE model
    %
    % offsetFree = false: 6 states Xa = [X; D], the linearized SDE (5.22)-(5.27)
    %   dX = (A*X + B*U + E*D)*dt,   dD = -Lambda*D*dt + Sigma*dw
    % offsetFree = true: 8 states [X; D; N], the mean-shift model (6.12)-(6.13)
    %   dD = -Lambda*(D - N)*dt + Sigma*dw,   dN = SigmaN*dwN
    % The model is discretized with c2dsde (ZOH + Van Loan), as in Chapter 5.
    %
    % Inputs:
    %   aug: Augmented SDE model created by computeAugmentedModel()
    %   lin: Linearized model created by computeJacobians()
    %   p: Parameter structure created by createParameters()
    %   Ts: Sampling time [s]
    %   offsetFree: true for the offset-free (8-state) model
    %   sigmaN: Diffusion of the mean shift N [cm^3/s/sqrt(s)] (tuning, Q_N)
    %
    % Output:
    %   model: Structure with the fields A, B, G, Q, C, R, S (discrete time, S = 0),
    %       P0, stateNames, stateUnits, scale (state -> physical units), name

    if nargin < 5, offsetFree = false; end
    if nargin < 6, sigmaN = [0.3; 0.3]; end

    A = aug.A;
    B = aug.B;
    G = aug.G;
    C = aug.C;

    stateNames = {'h_1', 'h_2', 'h_3', 'h_4', 'F_3', 'F_4'};
    stateUnits = {'cm', 'cm', 'cm', 'cm', 'cm^3/s', 'cm^3/s'};
    scale = [1 ./ (p.rho * p.A); 1; 1]; % masses -> levels, flows unchanged
    P0 = [(5 * p.rho * p.A) .^ 2; 20^2; 20^2];

    if offsetFree
        Lambda = diag(p.disturbanceRate);
        Gamma = [zeros(4, 2); Lambda];          % pulls D towards N
        A = [A, Gamma; zeros(2, 8)];            % N is a Brownian motion
        B = [B; zeros(2, 2)];
        G = [G, zeros(6, 2); zeros(2, 2), diag(sigmaN)];
        C = [C, zeros(2, 2)];
        stateNames = [stateNames, {'N_3', 'N_4'}];
        stateUnits = [stateUnits, {'cm^3/s', 'cm^3/s'}];
        scale = [scale; 1; 1];
        P0 = [P0; 20^2; 20^2];
    end

    [model.A, model.B, model.Q] = c2dsde(A, B, G, Ts);
    model.Q = (model.Q + model.Q') / 2;
    model.G = eye(size(model.A, 1));
    model.C = C;
    model.R = p.measurementCovariance;
    model.S = zeros(size(model.G, 2), size(C, 1)); % w and v independent
    model.P0 = diag(P0);
    model.stateNames = stateNames;
    model.stateUnits = stateUnits;
    model.scale = scale;
    model.offsetFree = offsetFree;
    model.problem = 5;
    model.lin = lin;
    model.name = sprintf('Problem 5 model (%d states)', size(model.A, 1));
end
