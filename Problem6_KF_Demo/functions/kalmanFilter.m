%% 
function result = kalmanFilter(model, y, u, filterType, x0, P0)
    if nargin < 4
        fprintf('Usage: result = kalmanFilter(model, y, u, ''dynamic'' or ''static''). Run demoKalmanFilter for an example.\n');
        if nargout > 0, result = []; end
        return
    end
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
