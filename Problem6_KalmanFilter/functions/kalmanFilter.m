function result = kalmanFilter(model, y, u, filterType, x0, P0)
    % Generic discrete-time Kalman filter, static or dynamic
    %
    % Works for ANY linear model given as a structure (lecture 6, slides 19-23):
    %   x(k+1) = A*x(k) + B*u(k) + G*w(k)
    %   y(k)   = C*x(k) + v(k)
    %   [w; v] ~ N(0, [Q S; S' R])      (S = cross-covariance, default 0)
    %
    % Dynamic filter (time-varying gain):
    %   Innovation:          e = y - C*xp,           Re = C*P*C' + R
    %   Gains:               Kfx = P*C'/Re,          Kfw = S/Re
    %   Measurement update:  xf = xp + Kfx*e,        wf = Kfw*e
    %                        Pf = P - Kfx*Re*Kfx',   Qf = Q - Kfw*Re*Kfw'
    %   Time update:         xp+ = A*xf + B*u + G*wf
    %                        P+  = A*Pf*A' + G*Qf*G' - A*Kfx*S'*G' - G*S*Kfx'*A'
    % Static filter: the same equations with the constant gains from
    % kalmanStatic() (solution of the discrete algebraic Riccati equation).
    %
    % Inputs:
    %   model: Structure with the fields
    %       A, C, Q, R            (required)
    %       B                     (optional, default: no input)
    %       G                     (optional, default: identity)
    %       S                     (optional, default: zeros)
    %       P0                    (optional, default: identity)
    %   y: ny x N measurements
    %   u: nu x N inputs ([] if the model has no input)
    %   filterType: 'dynamic' or 'static'
    %   x0: x_hat(0|-1) (optional, default: zeros)
    %   P0: P(0|-1)     (optional, default: model.P0)
    %
    % Output:
    %   result: Structure with the fields
    %       xf:  Filtered states x_hat(k|k), n x N
    %       xp:  One-step predictions x_hat(k|k-1), n x N
    %       wf:  Filtered process noise w_hat(k|k), nw x N
    %       e:   Innovations, ny x N
    %       Re:  Innovation covariances, ny x ny x N
    %       K:   Filter gains Kfx, n x ny x N
    %       sdf: Standard deviations of x_hat(k|k), sqrt(diag(Pf)), n x N
    %       sdp: Standard deviations of x_hat(k|k-1), sqrt(diag(P)), n x N

    model = completeModel(model);
    n = size(model.A, 1);
    ny = size(model.C, 1);
    nw = size(model.G, 2);
    N = size(y, 2);
    if isempty(u)
        u = zeros(size(model.B, 2), N);
    end
    if nargin < 5 || isempty(x0), x0 = zeros(n, 1); end
    if nargin < 6 || isempty(P0), P0 = model.P0; end

    A = model.A;  B = model.B;  G = model.G;  C = model.C;
    Q = model.Q;  R = model.R;  S = model.S;

    isStatic = strcmpi(filterType, 'static');
    if ~isStatic && ~strcmpi(filterType, 'dynamic')
        error('kalmanFilter:type', 'filterType must be ''static'' or ''dynamic''.');
    end
    if isStatic
        [~, KfxS, ReS, PfS, ~, KfwS, QfS] = kalmanStatic(model);
    end

    result.xf = zeros(n, N);
    result.xp = zeros(n, N);
    result.wf = zeros(nw, N);
    result.e = zeros(ny, N);
    result.Re = zeros(ny, ny, N);
    result.K = zeros(n, ny, N);
    result.sdf = zeros(n, N);
    result.sdp = zeros(n, N);

    xp = x0;
    P = P0;
    for k = 1:N
        % Gains
        if isStatic
            Re = ReS;  Kfx = KfxS;  Kfw = KfwS;  Pf = PfS;  Qf = QfS;
        else
            Re = C * P * C' + R;
            Re = (Re + Re') / 2;
            Kfx = P * C' / Re;
            Kfw = S / Re;
            Pf = P - Kfx * Re * Kfx';
            Pf = (Pf + Pf') / 2;
            Qf = Q - Kfw * Re * Kfw';
        end

        % Measurement update
        e = y(:, k) - C * xp;
        xf = xp + Kfx * e;
        wf = Kfw * e;

        result.xp(:, k) = xp;
        result.xf(:, k) = xf;
        result.wf(:, k) = wf;
        result.e(:, k) = e;
        result.Re(:, :, k) = Re;
        result.K(:, :, k) = Kfx;
        result.sdf(:, k) = sqrt(max(diag(Pf), 0));
        result.sdp(:, k) = sqrt(max(diag(P), 0));

        % Time update (one-step prediction)
        xp = A * xf + B * u(:, k) + G * wf;
        if ~isStatic
            P = A * Pf * A' + G * Qf * G' - A * Kfx * S' * G' - G * S * Kfx' * A';
            P = (P + P') / 2;
        end
    end

    result.filterType = lower(char(filterType));
    if isfield(model, 'name')
        result.modelName = char(model.name);
    else
        result.modelName = '';
    end
end

function model = completeModel(model)
    % Default values for the optional fields
    n = size(model.A, 1);
    ny = size(model.C, 1);
    if ~isfield(model, 'B') || isempty(model.B), model.B = zeros(n, 0); end
    if ~isfield(model, 'G') || isempty(model.G), model.G = eye(n); end
    if ~isfield(model, 'S') || isempty(model.S), model.S = zeros(size(model.G, 2), ny); end
    if ~isfield(model, 'P0') || isempty(model.P0), model.P0 = eye(n); end
end
