%% Problem 3 - PID control of the Modified Four Tank System
% Decentralised P, PI and PID control of the lower-tank levels h1 and h2.
% The plant is the stochastic nonlinear model with piecewise-constant
% disturbances F3, F4 and Gaussian measurement noise (Section 2.2).

clear; clc; close all;

%% Physical parameters
p.a     = 1.2272*ones(4,1);         % outlet cross-sectional areas [cm^2]
p.A     = 380.1327*ones(4,1);       % tank cross-sectional areas [cm^2]
p.gamma = [0.6; 0.7];               % valve split ratios [-]
p.g     = 981;                      % gravitational acceleration [cm/s^2]
p.rho   = 1.0;                      % density of water [g/cm^3]

%% Operating point (Table 5.1)
us = [300; 300];                    % nominal pump flows [cm^3/s]
ds = [50; 100];                     % nominal disturbance flows [cm^3/s]
[xs, hs, qs] = steadyState(us, ds, p);

%% Pairing analysis
[A, B, Cz] = linearise(hs, p);
Tc  = 2*p.A.*hs./qs;                % tank time constants [s]
G0  = -Cz*(A\B);                    % steady-state gain [cm/(cm^3/s)]
RGA = G0.*inv(G0).';
NI  = det(G0)/prod(diag(G0));       % Niederlinski index
z0  = transmissionZeros(A, B, Cz);

fprintf('Steady state h_s [cm]      : %s\n', mat2str(hs.', 4));
fprintf('Steady state q_s [cm^3/s]  : %s\n', mat2str(qs.', 4));
fprintf('Time constants T [s]       : %s\n', mat2str(Tc.', 4));
fprintf('Steady-state gain G0:\n'); disp(G0);
fprintf('Relative gain array:\n');  disp(RGA);
fprintf('lambda11 analytic = %.4f\n', prod(p.gamma)/(sum(p.gamma) - 1));
fprintf('Niederlinski index = %.4f\n', NI);
fprintf('Transmission zeros [1/s]   : %s\n\n', mat2str(z0.', 4));

%% Controller tuning
% The paired transfer functions g11 and g22 are first order,
% g_ii(s) = k_i/(T_i s + 1). The zero-order hold is approximated by a
% dead time theta = Ts/2. P and PI: SIMC rules. PID: IMC design on the
% first-order Pade approximation of the dead time, same closed-loop time constant.
Ts    = 10;                         % sampling time [s]
theta = Ts/2;
k     = diag(G0);
tau   = Tc(1:2);

tauC  = 4*theta;                    % desired closed-loop time constant
KcPI  = tau./(k*(tauC + theta));
TiPI  = min(tau, 4*(tauC + theta));

KcPID = (tau + theta/2)./(k.*(tauC + theta/2));
TiPID = tau + theta/2;
TdPID = tau*theta./(2*tau + theta);
TfPID = TdPID/10;                   % derivative filter, N = 10

ctrl(1) = struct('name','P',   'Kc',KcPI,  'Ti',[Inf;Inf], 'Td',[0;0],  'Tf',[0;0]);
ctrl(2) = struct('name','PI',  'Kc',KcPI,  'Ti',TiPI,      'Td',[0;0],  'Tf',[0;0]);
ctrl(3) = struct('name','PID', 'Kc',KcPID, 'Ti',TiPID,     'Td',TdPID,  'Tf',TfPID);

fprintf('%-4s | %8s %8s | %7s %7s | %6s %6s | %6s %6s\n', 'Ctrl', ...
        'Kc1','Kc2','Ti1','Ti2','Td1','Td2','Tf1','Tf2');
for c = 1:numel(ctrl)
    fprintf('%-4s | %8.2f %8.2f | %7.1f %7.1f | %6.2f %6.2f | %6.2f %6.2f\n', ...
        ctrl(c).name, ctrl(c).Kc, ctrl(c).Ti, ctrl(c).Td, ctrl(c).Tf);
end
fprintf('\n');

%% Stability of the linearised sampled closed loop (both loops closed)
Md = expm([A B; zeros(2,6)]*Ts);
Ad = Md(1:4,1:4);
Bd = Md(1:4,5:6);
ctrlChk = [ctrl, ctrl(3)];
ctrlChk(4).name = 'PID, Td = 0';
ctrlChk(4).Td = [0; 0];
fprintf('%-12s | %8s | %10s\n', 'Ctrl', 'rho(Acl)', 'GM factor');
for c = 1:numel(ctrlChk)
    rhoCl = max(abs(eig(closedLoopMatrix(Ad, Bd, Cz, ctrlChk(c), Ts))));
    gm = NaN;
    cs = ctrlChk(c);
    for f = 1:0.01:20
        cs.Kc = f*ctrlChk(c).Kc;
        if max(abs(eig(closedLoopMatrix(Ad, Bd, Cz, cs, Ts)))) >= 1
            gm = f;
            break
        end
    end
    fprintf('%-12s | %8.4f | %10.2f\n', ctrlChk(c).name, rhoCl, gm);
end
fprintf('(GM factor: common scaling of both gains Kc at which the loop becomes unstable)\n\n');

%% Predicted steady-state offsets of the P-controller (linear model)
E   = p.rho*[0 0; 0 0; 1 0; 0 1];
Gd0 = -Cz*(A\E);
S0  = inv(eye(2) + G0*diag(KcPI));
fprintf('Linear P offsets e = S0*(dr - Gd0*dd) [cm]\n');
fprintf('  after r1 step         : %s\n', mat2str((S0*[10; 0]).', 3));
fprintf('  after r1 and r2 steps : %s\n', mat2str((S0*[10; -10]).', 3));
fprintf('  after F3 step         : %s\n\n', mat2str((S0*([10; -10] - Gd0*[30; 0])).', 3));

%% Test scenario
tEnd  = 40*60;
t     = (0:Ts:tEnd)';
N     = numel(t);
uMin  = [0; 0];
uMax  = [500; 500];                 % pump capacity [cm^3/s]

r = repmat(hs(1:2), 1, N);
r(1, t >= 5*60)  = hs(1) + 10;      % setpoint step in loop 1
r(2, t >= 15*60) = hs(2) - 10;      % setpoint step in loop 2

dBar = repmat(ds, 1, N);
dBar(1, t >= 25*60) = ds(1) + 30;   % unmeasured step in F3

rng(100, 'twister');
sigD = [10; 20];                    % std. dev. of F3, F4 [cm^3/s]
Rvv  = eye(2);                      % measurement noise covariance [cm^2]
d    = max(0, dBar + sigD.*randn(2, N));
v    = chol(Rvv, 'lower')*randn(2, N);

%% Closed-loop simulations
% Deterministic case: mean disturbances, no noise. Stochastic case: model (2.2).
for c = 1:numel(ctrl)
    resDet(c) = simulateClosedLoop(ctrl(c), xs, us, t, r, dBar, 0*v, uMin, uMax, p); %#ok<SAGROW>
    resSto(c) = simulateClosedLoop(ctrl(c), xs, us, t, r, d, v, uMin, uMax, p);       %#ok<SAGROW>
end

%% Performance metrics
fprintf('Deterministic case\n');
fprintf('%-4s | %6s %6s | %6s %6s | %6s %6s | %6s %6s | %6s\n', 'Ctrl', 'IAE1', 'IAE2', ...
        't90', 'ts', 'e1', 'e2', 'eD1', 'eD2', 'dh2');
for c = 1:numel(resDet)
    resDet(c).metrics = trackingMetrics(t, r, resDet(c));
    fprintf('%-4s | %6.1f %6.1f | %6.2f %6.2f | %6.2f %6.2f | %6.2f %6.2f | %6.2f\n', ...
            resDet(c).name, resDet(c).metrics);
end
fprintf('(IAE [cm min]; t90, ts [min] after the h1 step; e1, e2, eD1, eD2 [cm]; dh2 [cm])\n\n');

fprintf('Stochastic case\n');
fprintf('%-4s | %6s %6s | %6s %6s | %6s %6s\n', 'Ctrl', 'IAE1', 'IAE2', ...
        'sd_e1', 'sd_e2', 'TV1', 'TV2');
iSS = t < 5*60;
for c = 1:numel(resSto)
    e   = r - resSto(c).h(1:2,:);
    IAE = sum(abs(e), 2)*Ts/60;
    sdE = std(e(:, iSS), 0, 2);
    TV  = mean(abs(diff(resSto(c).u, 1, 2)), 2);
    resSto(c).metrics = [IAE.' sdE.' TV.'];
    fprintf('%-4s | %6.1f %6.1f | %6.2f %6.2f | %6.1f %6.1f\n', resSto(c).name, resSto(c).metrics);
end
fprintf('(IAE [cm min]; sd_e: std of e over 0-5 min [cm]; TV = mean|du_k| [cm^3/s])\n\n');

%% Figures
col = [0.0000 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.4660 0.6740 0.1880];
outDir = 'figures';
if ~exist(outDir, 'dir'), mkdir(outDir); end

plotDisturbances(t, dBar, d, fullfile(outDir, 'P3_disturbances.png'));
for c = 1:numel(resSto)
    plotLoops(t, r, resSto(c), col(c,:), ...
        sprintf('Problem 3: %s-Controller - Closed Loop', resSto(c).name), ...
        fullfile(outDir, sprintf('P3_%s.png', resSto(c).name)));
end
plotComparison(t, r, resDet, col, [0 40], ...
    'Problem 3: P, PI and PID - Deterministic Model', ...
    fullfile(outDir, 'P3_comparison_deterministic.png'));
plotComparison(t, r, resDet, col, [4 20], ...
    'Problem 3: P, PI and PID - Setpoint Steps (Deterministic)', ...
    fullfile(outDir, 'P3_comparison_zoom.png'));
plotComparison(t, r, resSto, col, [0 40], ...
    'Problem 3: P, PI and PID - Stochastic Model', ...
    fullfile(outDir, 'P3_comparison_stochastic.png'));
plotUpperTanks(t, resDet, col, fullfile(outDir, 'P3_upper_tanks.png'));

%% ------------------------------------------------------------------------
%  Local functions
%  ------------------------------------------------------------------------
function s = simulateClosedLoop(c, x0, ubar, t, r, d, v, uMin, uMax, p)
    % Discrete-time decentralised control of the continuous nonlinear plant.
    % Inputs and disturbances are held constant over each sample (ZOH).
    N  = numel(t);
    Ts = t(2) - t(1);
    x  = x0;
    h  = zeros(4, N);
    y  = zeros(2, N);
    u  = zeros(2, N);
    st = struct('I', [0;0], 'D', [0;0], 'yPrev', []);
    opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);
    for kk = 1:N
        h(:,kk) = x./(p.rho*p.A);
        y(:,kk) = h(1:2,kk) + v(:,kk);
        [u(:,kk), st] = pidControl(r(:,kk), y(:,kk), st, c, ubar, uMin, uMax, Ts);
        if kk < N
            [~, X] = ode45(@(tt, xx) modifiedFourTank(xx, u(:,kk), d(:,kk), p), ...
                           [t(kk) t(kk+1)], x, opts);
            x = max(X(end,:).', 0);
        end
    end
    s = struct('name', c.name, 'h', h, 'y', y, 'u', u, 'metrics', []);
end

function m = trackingMetrics(t, r, s)
    % Tracking and regulation metrics of a deterministic closed-loop run.
    Ts  = t(2) - t(1);
    tm  = t/60;
    e   = r - s.h(1:2,:);
    IAE = sum(abs(e), 2)*Ts/60;
    iR1 = tm >= 5 & tm < 15;
    h1  = s.h(1, iR1);
    t1  = tm(iR1) - 5;
    dr  = r(1, find(iR1, 1)) - r(1, find(iR1, 1) - 1);
    h10 = s.h(1, find(iR1, 1) - 1);
    i90 = find(h1 - h10 >= 0.9*dr, 1);
    out = find(abs(r(1, iR1) - h1) > 0.05*abs(dr), 1, 'last');
    t90 = NaN; ts = NaN;
    if ~isempty(i90), t90 = t1(i90); end
    if isempty(out), ts = 0; elseif out < numel(t1), ts = t1(out + 1); end
    if abs(e(1, find(iR1, 1, 'last'))) > 0.05*abs(dr), ts = NaN; end
    e1  = e(1, find(tm < 15, 1, 'last'));
    e2  = e(2, find(tm < 25, 1, 'last'));
    eD  = e(:, end);
    dh2 = max(abs(e(2, iR1)));
    m   = [IAE.' t90 ts e1 e2 eD.' dh2];
end

function Acl = closedLoopMatrix(Ad, Bd, Cz, c, Ts)
    % Linearised sampled closed loop, unsaturated controller (3.10a)-(3.10e).
    % States: [x_k; I_k; D_{k-1}; y_{k-1}] in deviation variables.
    hasI = all(isfinite(c.Ti));
    hasD = any(c.Td > 0);
    Kc = diag(c.Kc);
    ki = zeros(2); a = zeros(2); b = zeros(2);
    if hasI, ki = diag(c.Kc*Ts./c.Ti); end
    if hasD
        a = diag(c.Tf./(c.Tf + Ts));
        b = diag(c.Kc.*c.Td./(c.Tf + Ts));
    end
    Acl = [Ad - Bd*(Kc + b)*Cz, Bd,       Bd*a,     Bd*b;
           -ki*Cz,              eye(2),   zeros(2), zeros(2);
           -b*Cz,               zeros(2), a,        b;
           Cz,                  zeros(2), zeros(2), zeros(2)];
    keep = [true(4,1); repmat(hasI, 2, 1); repmat(hasD, 4, 1)];
    Acl = Acl(keep, keep);
end

function [u, st] = pidControl(r, y, st, c, ubar, uMin, uMax, Ts)
    % Positional PID with derivative on the measurement,
    % first-order derivative filter and conditional-integration anti-windup.
    e = r - y;
    if any(c.Td > 0) && ~isempty(st.yPrev)
        st.D = c.Tf./(c.Tf + Ts).*st.D - c.Kc.*c.Td./(c.Tf + Ts).*(y - st.yPrev);
    end
    st.yPrev = y;
    v = ubar + c.Kc.*e + st.I + st.D;
    u = min(max(v, uMin), uMax);
    if all(isfinite(c.Ti))
        integrate = (v == u) | (sign(e) ~= sign(v - u));
        st.I = st.I + integrate.*c.Kc*Ts./c.Ti.*e;
    end
end

function xdot = modifiedFourTank(x, u, d, p)
    % Mass balances of the Modified Four Tank System, eq. (2.10).
    h    = max(x, 0)./(p.rho*p.A);
    q    = p.a.*sqrt(2*p.g*h);
    g1   = p.gamma(1);
    g2   = p.gamma(2);
    qin  = [g1*u(1); g2*u(2); (1-g2)*u(2); (1-g1)*u(1)];
    xdot = p.rho*[qin(1) + q(3) - q(1);
                  qin(2) + q(4) - q(2);
                  qin(3) - q(3) + d(1);
                  qin(4) - q(4) + d(2)];
end

function [xs, hs, qs] = steadyState(us, ds, p)
    % Analytical steady state from the flow balances, eq. (5.6).
    g1 = p.gamma(1);
    g2 = p.gamma(2);
    q3 = (1-g2)*us(2) + ds(1);
    q4 = (1-g1)*us(1) + ds(2);
    qs = [g1*us(1) + q3; g2*us(2) + q4; q3; q4];
    hs = qs.^2./(2*p.g*p.a.^2);
    xs = p.rho*p.A.*hs;
end

function [A, B, Cz] = linearise(hs, p)
    % Jacobians of the drift at the steady state, eqs. (5.10)-(5.11).
    T  = 2*p.A.*hs./(p.a.*sqrt(2*p.g*hs));
    A  = diag(-1./T);
    A(1,3) = 1/T(3);
    A(2,4) = 1/T(4);
    g1 = p.gamma(1);
    g2 = p.gamma(2);
    B  = p.rho*[g1 0; 0 g2; 0 1-g2; 1-g1 0];
    Cz = [1/(p.rho*p.A(1)) 0 0 0; 0 1/(p.rho*p.A(2)) 0 0];
end

function z = transmissionZeros(A, B, C)
    % Finite generalised eigenvalues of the Rosenbrock system matrix.
    n = size(A,1);
    m = size(B,2);
    M = [A B; -C zeros(m)];
    E = blkdiag(eye(n), zeros(m));
    z = eig(M, E);
    z = z(isfinite(z));
end

function plotLoops(t, r, s, c, ttl, file)
    fig = figure('Position', [100 100 900 750]);
    tm  = t/60;
    for j = 1:2
        subplot(3,2,j); hold on; grid on; box on;
        plot(tm, s.y(j,:).', '-', 'Color', 0.55*c + 0.45, 'LineWidth', 0.7);
        plot(tm, s.h(j,:).', '-', 'Color', c, 'LineWidth', 1.6);
        stairs(tm, r(j,:).', 'k--', 'LineWidth', 1.2);
        xlabel('Time [min]'); ylabel(sprintf('h_%d [cm]', j));
        title(sprintf('Tank %d level', j));
        if j == 1
            legend({'measured', 'true', 'setpoint'}, 'Location', 'southeast');
        end

        subplot(3,2,2+j); hold on; grid on; box on;
        stairs(tm, s.u(j,:).', '-', 'Color', c, 'LineWidth', 1.2);
        xlabel('Time [min]'); ylabel(sprintf('F_%d [cm^3/s]', j));
        title(sprintf('Pump %d flow', j)); ylim([-20 520]);

        subplot(3,2,4+j); hold on; grid on; box on;
        plot(tm, (r(j,:) - s.h(j,:)).', '-', 'Color', c, 'LineWidth', 1.4);
        plot(tm([1 end]), [0 0], 'k:');
        xlabel('Time [min]'); ylabel(sprintf('e_%d [cm]', j));
        title(sprintf('Tracking error, tank %d', j));
    end
    sgtitle(ttl);
    print(fig, file, '-dpng', '-r200');
end

function plotComparison(t, r, res, col, win, ttl, file)
    fig = figure('Position', [100 100 900 750]);
    tm  = t/60;
    for j = 1:2
        subplot(3,2,j); hold on; grid on; box on;
        for c = 1:numel(res)
            plot(tm, res(c).h(j,:).', 'Color', col(c,:), 'LineWidth', 1.5);
        end
        stairs(tm, r(j,:).', 'k--', 'LineWidth', 1.2);
        xlim(win); xlabel('Time [min]'); ylabel(sprintf('h_%d [cm]', j));
        title(sprintf('Tank %d level', j));
        if j == 1
            legend([{res.name}, {'setpoint'}], 'Location', 'southeast');
        end

        subplot(3,2,2+j); hold on; grid on; box on;
        for c = 1:numel(res)
            stairs(tm, res(c).u(j,:).', 'Color', col(c,:), 'LineWidth', 1.1);
        end
        xlim(win); ylim([-20 520]);
        xlabel('Time [min]'); ylabel(sprintf('F_%d [cm^3/s]', j));
        title(sprintf('Pump %d flow', j));

        subplot(3,2,4+j); hold on; grid on; box on;
        for c = 1:numel(res)
            plot(tm, (r(j,:) - res(c).h(j,:)).', 'Color', col(c,:), 'LineWidth', 1.4);
        end
        plot(win, [0 0], 'k:');
        xlim(win); xlabel('Time [min]'); ylabel(sprintf('e_%d [cm]', j));
        title(sprintf('Tracking error, tank %d', j));
    end
    sgtitle(ttl);
    print(fig, file, '-dpng', '-r200');
end

function plotDisturbances(t, dBar, d, file)
    fig = figure('Position', [100 100 900 450]);
    tm  = t/60;
    c   = [0.0000 0.4470 0.7410; 0.8500 0.3250 0.0980];
    for j = 1:2
        subplot(2,1,j); hold on; grid on; box on;
        stairs(tm, d(j,:).', 'Color', c(j,:), 'LineWidth', 1.1);
        stairs(tm, dBar(j,:).', 'k--', 'LineWidth', 1.3);
        xlabel('Time [min]'); ylabel(sprintf('F_%d [cm^3/s]', j+2));
        if j == 1, legend({'realisation d_k', 'mean'}, 'Location', 'northwest'); end
    end
    sgtitle('Problem 3: Piecewise-Constant Disturbance Flows');
    print(fig, file, '-dpng', '-r200');
end

function plotUpperTanks(t, res, col, file)
    fig = figure('Position', [100 100 900 380]);
    tm  = t/60;
    for j = 1:2
        subplot(1,2,j); hold on; grid on; box on;
        for c = 1:numel(res)
            plot(tm, res(c).h(j+2,:).', 'Color', col(c,:), 'LineWidth', 1.4);
        end
        xlabel('Time [min]'); ylabel(sprintf('h_%d [cm]', j+2));
        title(sprintf('Tank %d level', j+2));
        if j == 1, legend({res.name}, 'Location', 'northwest'); end
    end
    sgtitle('Problem 3: Upper Tank Levels Under P, PI and PID Control');
    print(fig, file, '-dpng', '-r200');
end
