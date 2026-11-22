function model = buildKalmanModelP4(ident, p, Ts, offsetFree, sigmaW, sigmaN)
    % Discrete-time Kalman filter model from the identified model of Problem 4
    %
    % Minimal 4-state realization of (4.12) with shared T1 and T2:
    %   x1' = (-x1 + K11*u1 + K12*x3)/T1      y1 = x1
    %   x2' = (-x2 + K22*u2 + K21*x4)/T2      y2 = x2
    %   x3' = (-x3 + u2 + w3)/T3
    %   x4' = (-x4 + u1 + w4)/T4
    % x1, x2 are level deviations [cm]. x3, x4 are filtered flows in pump-flow
    % units [cm^3/s] and play the role of the upper tanks.
    %
    % The identified model has no noise description. The noise is ASSUMED to
    % enter at x3 and x4, where F3 and F4 enter physically, and to be held
    % constant over a sample like the inputs: w ~ N(0, diag(sigmaW.^2)).
    % offsetFree = true adds integrating disturbances N on the same channel:
    %   x(k+1) = Ad*x + Bd*u + Bwd*(N + w),   N(k+1) = N(k) + xi,
    %   xi ~ N(0, diag(sigmaN.^2))  -> 6 states.
    %
    % Inputs:
    %   ident: Structure with K (2x2) and T (4x1), see sharedTimeConstants()
    %   p: Parameter structure created by createParameters()
    %   Ts: Sampling time [s]
    %   offsetFree: true for the offset-free (6-state) model
    %   sigmaW, sigmaN: Tuning of the noise and of the integrators [cm^3/s]
    %
    % Output:
    %   model: Structure with the fields A, B, G, Q, C, R, S (discrete time, S = 0),
    %       P0, stateNames, stateUnits, scale, name

    if nargin < 4, offsetFree = false; end
    if nargin < 5, sigmaW = [30; 30]; end
    if nargin < 6, sigmaN = [3; 3]; end

    K = ident.K;
    T = ident.T;

    A = [-1/T(1), 0,       K(1, 2)/T(1), 0;
          0,      -1/T(2), 0,            K(2, 1)/T(2);
          0,       0,      -1/T(3),      0;
          0,       0,       0,           -1/T(4)];
    B = [K(1, 1)/T(1), 0;
         0,            K(2, 2)/T(2);
         0,            1/T(3);
         1/T(4),       0];
    Bw = [0, 0; 0, 0; 1/T(3), 0; 0, 1/T(4)];
    C = [1, 0, 0, 0; 0, 1, 0, 0];

    % ZOH discretization of the inputs and the noise together
    [Ad, BBw] = c2dzoh(A, [B, Bw], Ts);
    Bd = BBw(:, 1:2);
    Bwd = BBw(:, 3:4);

    stateNames = {'x_1 (h_1)', 'x_2 (h_2)', 'x_3', 'x_4'};
    stateUnits = {'cm', 'cm', 'cm^3/s', 'cm^3/s'};

    if ~offsetFree
        model.A = Ad;
        model.B = Bd;
        model.G = Bwd;
        model.Q = diag(sigmaW .^ 2);
        model.C = C;
        P0 = [25; 25; 50^2; 50^2];
    else
        model.A = [Ad, Bwd; zeros(2, 4), eye(2)];
        model.B = [Bd; zeros(2, 2)];
        model.G = [Bwd, zeros(4, 2); zeros(2, 2), eye(2)];
        model.Q = blkdiag(diag(sigmaW .^ 2), diag(sigmaN .^ 2));
        model.C = [C, zeros(2, 2)];
        stateNames = [stateNames, {'N_1', 'N_2'}];
        stateUnits = [stateUnits, {'cm^3/s', 'cm^3/s'}];
        P0 = [25; 25; 50^2; 50^2; 50^2; 50^2];
    end

    model.R = p.measurementCovariance;
    model.S = zeros(size(model.G, 2), size(model.C, 1)); % w and v independent
    model.P0 = diag(P0);
    model.stateNames = stateNames;
    model.stateUnits = stateUnits;
    model.scale = ones(size(model.A, 1), 1);
    model.offsetFree = offsetFree;
    model.problem = 4;
    model.ident = ident;
    model.name = sprintf('Problem 4 model (%d states)', size(model.A, 1));
end
