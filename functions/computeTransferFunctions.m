function trf = computeTransferFunctions(lin, p)
    % Continuous-time transfer functions of the linearized four-tank system
    %
    %   Z(s) = G(s)*U(s) + Gd(s)*D(s)              (deterministic, piecewise)
    %   Z(s) = G(s)*U(s) + H(s)*dw(s)              (SDE)
    % with G(s) = Cz*(sI - A)^(-1)*B, Gd(s) = Cz*(sI - A)^(-1)*E and
    % H(s) = Gd(s)*diag(sigma_F./(s + a_F)).
    %
    % Every entry is a gain times first-order lags: K_ij / prod_l (T_l*s + 1).
    %
    % Inputs:
    %   lin: Linearized model created by computeJacobians()
    %   p: Parameter structure created by createParameters()
    %
    % Output:
    %   trf: Structure with the fields
    %       G, Gd, H: Function handles s -> 2x2 transfer function matrix
    %       Ku, Tu: Gains [cm/(cm^3/s)] and time constants [s] of G(s)
    %       Kd, Td: Gains and time constants of Gd(s)
    %       Kh, Th: Gains and time constants of H(s)
    %   Time constants are stored in 2x2 cells; an empty entry means a zero entry.

    gamma = p.gamma;            % Valve positions [-]
    At = p.A;                   % Tank areas [cm^2]
    T = lin.T;                  % Tank time constants [s]
    a = p.disturbanceRate;      % Mean-reversion rates [1/s]
    sigma = p.disturbanceSigma; % Diffusion coefficients

    % Transfer function matrices, evaluated at a given (complex) s
    I = eye(size(lin.A));
    trf.G = @(s) lin.Cz / (s * I - lin.A) * lin.B;
    trf.Gd = @(s) lin.Cz / (s * I - lin.A) * lin.E;
    trf.H = @(s) trf.Gd(s) * diag(sigma ./ (s + a));

    % Steady-state level change of tanks 1 and 2 per unit inflow [cm/(cm^3/s)]
    k = [T(1) / At(1); T(2) / At(2)];

    % From the pumps: direct path (one lag) or via the top tank (two lags)
    trf.Ku = [gamma(1) * k(1),       (1 - gamma(2)) * k(1);
              (1 - gamma(1)) * k(2), gamma(2) * k(2)];
    trf.Tu = {T(1),          [T(1), T(3)];
              [T(2), T(4)],  T(2)};

    % From the disturbances: always via the top tanks
    trf.Kd = diag(k);
    trf.Td = {[T(1), T(3)], [];
              [],           [T(2), T(4)]};

    % From the Brownian noise: the disturbance filter adds the lag 1/a
    trf.Kh = trf.Kd * diag(sigma ./ a);
    trf.Th = {[trf.Td{1, 1}, 1 / a(1)], [];
              [],                    [trf.Td{2, 2}, 1 / a(2)]};
end
