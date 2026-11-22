function [P, Kfx, Re, Pf, info, Kfw, Qf] = kalmanStatic(model)
    % Stationary (static) Kalman filter for any linear model (lecture 6, slides 28-30)
    %
    % Model: x(k+1) = A*x(k) + B*u(k) + G*w(k),  y(k) = C*x(k) + v(k),
    %        [w; v] ~ N(0, [Q S; S' R])  (S optional, default 0)
    %
    % Discrete algebraic Riccati equation for the prediction covariance P:
    %   P = A*P*A' + G*Q*G' - (A*P*C' + G*S)*(C*P*C' + R)^(-1)*(A*P*C' + G*S)'
    % Solved with idare/dare (slide 29, method 1) when the Control System
    % Toolbox is available, and always also by fixed-point iteration of the
    % Riccati recursion (slide 29, method 3) as a check. Without the toolbox
    % the iteration result is used.
    %
    % Outputs:
    %   P:   Stationary one-step prediction covariance P(k+1|k)
    %   Kfx: Filter gain, Kfx = P*C'/Re
    %   Re:  Innovation covariance, Re = C*P*C' + R
    %   Pf:  Filtered covariance, Pf = P - Kfx*Re*Kfx'
    %   info: method, iterations, relDiffIteration (solver vs iteration),
    %         eigErrorDynamics (eigenvalues of A - Kpx*C, all inside the unit circle)
    %   Kfw: Noise gain, Kfw = S/Re
    %   Qf:  Filtered noise covariance, Qf = Q - Kfw*Re*Kfw'

    A = model.A;
    C = model.C;
    Q = model.Q;
    R = model.R;
    n = size(A, 1);
    ny = size(C, 1);
    if isfield(model, 'G') && ~isempty(model.G), G = model.G; else, G = eye(n); end
    if isfield(model, 'S') && ~isempty(model.S), S = model.S; else, S = zeros(size(G, 2), ny); end
    GQG = G * Q * G';
    GS = G * S;

    % Method 3: fixed-point iteration of the Riccati recursion
    Pit = GQG;
    info.iterations = 0;
    for it = 1:500000
        M = A * Pit * C' + GS;
        Pnew = A * Pit * A' + GQG - M / (C * Pit * C' + R) * M';
        Pnew = (Pnew + Pnew') / 2;
        if norm(Pnew - Pit, 'fro') <= 1e-12 * max(1, norm(Pnew, 'fro'))
            Pit = Pnew;
            info.iterations = it;
            break
        end
        Pit = Pnew;
    end
    if info.iterations == 0
        warning('kalmanStatic:iteration', 'Riccati iteration did not converge.');
    end

    % Method 1: idare / dare (Control System Toolbox)
    P = [];
    info.method = 'fixed-point iteration';
    if exist('idare', 'file') == 2 || exist('idare', 'builtin') == 5
        try
            [P, ~, ~, report] = idare(A', C', GQG, R, GS);
            if isempty(P) || (isstruct(report) && isfield(report, 'Report') && report.Report ~= 0)
                P = [];
            else
                info.method = 'idare';
            end
        catch
            P = [];
        end
    end
    if isempty(P) && (exist('dare', 'file') == 2 || exist('dare', 'builtin') == 5)
        try
            P = dare(A', C', GQG, R, GS);
            info.method = 'dare';
        catch
            P = [];
        end
    end
    if isempty(P)
        P = Pit;
        info.method = 'fixed-point iteration';
    end
    info.relDiffIteration = norm(P - Pit, 'fro') / max(norm(P, 'fro'), eps);

    P = (P + P') / 2;
    Re = C * P * C' + R;
    Re = (Re + Re') / 2;
    Kfx = P * C' / Re;
    Kfw = S / Re;
    Pf = P - Kfx * Re * Kfx';
    Pf = (Pf + Pf') / 2;
    Qf = Q - Kfw * Re * Kfw';
    Kpx = A * Kfx + G * Kfw;                % one-step predictor gain
    info.eigErrorDynamics = eig(A - Kpx * C);
end
