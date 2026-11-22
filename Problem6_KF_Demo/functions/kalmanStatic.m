function [P, Kfx, Re, Pf, info, Kfw, Qf] = kalmanStatic(model)
% KALMANSTATIC  Stationary Kalman filter from the discrete algebraic Riccati equation.
    A = model.A;  C = model.C;  Q = model.Q;  R = model.R;
    n = size(A, 1);  ny = size(C, 1);
    if isfield(model, 'G') && ~isempty(model.G), G = model.G; else, G = eye(n); end
    if isfield(model, 'S') && ~isempty(model.S), S = model.S; else, S = zeros(size(G, 2), ny); end
    GQG = G * Q * G';
    GS = G * S;

    % Fixed-point iteration of the Riccati recursion
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
    Kpx = A * Kfx + G * Kfw;
    info.eigErrorDynamics = eig(A - Kpx * C);
end
