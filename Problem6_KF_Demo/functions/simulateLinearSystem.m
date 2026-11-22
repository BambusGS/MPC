function sim = simulateLinearSystem(model, u, N, x0mean)
% SIMULATELINEARSYSTEM  Simulate x(k+1) = A x + B u + G w, y = C x + v, [w; v] ~ N(0, [Q S; S' R]).
    n = size(model.A, 1);
    ny = size(model.C, 1);
    if ~isfield(model, 'B') || isempty(model.B), model.B = zeros(n, 0); end
    if ~isfield(model, 'G') || isempty(model.G), model.G = eye(n); end
    nw = size(model.G, 2);
    if ~isfield(model, 'S') || isempty(model.S), model.S = zeros(nw, ny); end
    if ~isfield(model, 'P0') || isempty(model.P0), model.P0 = eye(n); end
    if isempty(u), u = zeros(size(model.B, 2), N); end
    if nargin < 4 || isempty(x0mean), x0mean = zeros(n, 1); end

    Sigma = [model.Q, model.S; model.S', model.R];
    Sigma = (Sigma + Sigma') / 2;
    L = chol(Sigma + 1e-12 * eye(nw + ny), 'lower');

    sim.x = zeros(n, N);
    sim.y = zeros(ny, N);
    sim.w = zeros(nw, N);
    sim.v = zeros(ny, N);
    x = x0mean + chol(model.P0 + 1e-12 * eye(n), 'lower') * randn(n, 1);
    for k = 1:N
        wv = L * randn(nw + ny, 1);
        sim.w(:, k) = wv(1:nw);
        sim.v(:, k) = wv(nw + 1:end);
        sim.x(:, k) = x;
        sim.y(:, k) = model.C * x + sim.v(:, k);
        x = model.A * x + model.B * u(:, k) + model.G * sim.w(:, k);
    end
end
