% FOUR_TANK_CONTROL  Input/output pairing, P / PI / PID control and

clear; clc; close all;

%% ------------------------------------------------------------------------
%  PARAMETERS 
% -------------------------------------------------------------------------
a1 = 1.2272; a2 = 1.2272; a3 = 1.2272; a4 = 1.2272;          % [cm^2]
A1 = 380.1327; A2 = 380.1327; A3 = 380.1327; A4 = 380.1327;  % [cm^2]
gamma1 = 0.65; gamma2 = 0.55;                                % [-]
g = 981; rho = 1.00;                                         % [cm/s^2], [g/cm^3]
p = [a1; a2; a3; a4; A1; A2; A3; A4; gamma1; gamma2; g; rho];

%% ------------------------------------------------------------------------
%  1. PAIRING OF INPUTS AND OUTPUTS
% -------------------------------------------------------------------------

F_op = [250; 325];                    
h_op = steadyStateLevels(F_op, p);     % closed form, no iterative solver needed
fprintf('Steady-state levels h1..h4 [cm]: %s\n', mat2str(h_op.', 4));

[A_lin, B_lin, C_lin] = linearise(h_op, p);
C_z = C_lin(1:2, :);
G0  = -C_z * (A_lin \ B_lin);           % steady-state gain z = G0*u  [cm/(cm^3/s)]
RGA = G0 .* inv(G0).';
lambda_analytic = gamma1*gamma2 / (gamma1 + gamma2 - 1);

fprintf('--- Pairing analysis at F = [%g; %g] cm^3/s ---\n', F_op);
fprintf('Steady-state gain G0 [cm per cm^3/s]:\n'); disp(G0);
fprintf('RGA:\n'); disp(RGA);
fprintf('lambda11 analytic = g1*g2/(g1+g2-1) = %.4f\n', lambda_analytic);
fprintf('Linear zeros (F -> h1,h2): \n'); disp(tzero_2x2(A_lin, B_lin, C_z));
fprintf('gamma1+gamma2 = %.2f > 1  -> minimum-phase configuration\n', gamma1+gamma2);
if RGA(1,1) > 0.5
    fprintf('=> Pair F1 -> h1 and F2 -> h2 (diagonal pairing).\n\n');
else
    fprintf('=> Pair F1 -> h2 and F2 -> h1 (off-diagonal pairing).\n\n');
end

%% ------------------------------------------------------------------------
%  2. CONTROLLER SETTINGS
% -------------------------------------------------------------------------
Ts   = 10;            
t_f  = 40*60;              % simulation length [s]
t    = (0:Ts:t_f)';
N    = numel(t);

u_bias = [300; 300];       
u_min  = [0; 0];           % pumps cannot run backwards
u_max  = [500; 500];       % ASSUMED pump capacity [cm^3/s] - not given in source

% setpoint step on h1 at t = 20 min to expose interaction in loop 2.
r = repmat([30; 30], 1, N);
r(1, t >= 20*60) = 35;

% Tuning (same for both loops). Kc in [cm^3/s per cm], Ti/Td/Tf in [s].
ctrl(1) = struct('name', 'P',   'Kc', 20, 'Ti', Inf, 'Td', 0,  'Tf', 0);
ctrl(2) = struct('name', 'PI',  'Kc', 20, 'Ti', 150, 'Td', 0,  'Tf', 0);
ctrl(3) = struct('name', 'PID', 'Kc', 20, 'Ti', 150, 'Td', 10, 'Tf', 10);

%% ------------------------------------------------------------------------
%  NOISE 
% -------------------------------------------------------------------------
randn('state', 100);
Q = [20^2 0; 0 40^2];  w = chol(Q, 'lower') * randn(2, N);   % pump noise
R = eye(4);            v = chol(R, 'lower') * randn(4, N);   % level noise

%% ------------------------------------------------------------------------
%  3. CLOSED-LOOP SIMULATIONS
% -------------------------------------------------------------------------
x0 = zeros(4, 1);   
res = struct([]);
for c = 1:numel(ctrl)
    for useNoise = [false true]
        s = simulateClosedLoop(ctrl(c), p, x0, t, r, u_bias, u_min, u_max, ...
                               w*useNoise, v*useNoise);
        s.noise = useNoise;
        res = [res, s]; %#ok<AGROW>
    end
end

%% Performance table ------------------------------------------------------
fprintf('%-4s %-6s | %8s %8s | %8s %8s | %9s %9s | %8s\n', 'Ctrl', 'Noise', ...
        'e1(20m)', 'e2(20m)', 'e1(end)', 'e2(end)', 'IAE1', 'IAE2', 'max h2');
for k = 1:numel(res)
    s = res(k); i20 = find(t < 20*60, 1, 'last'); after = t >= 20*60;
    e = r - s.h(1:2, :);
    fprintf('%-4s %-6s | %8.2f %8.2f | %8.2f %8.2f | %9.1f %9.1f | %8.2f\n', ...
        s.name, mat2str(s.noise), e(1,i20), e(2,i20), e(1,end), e(2,end), ...
        sum(abs(e(1,:)))*Ts/60, sum(abs(e(2,:)))*Ts/60, max(s.h(2, after)));
end
fprintf('(e = r - h in cm, IAE in cm*min, max h2 after the h1 setpoint step)\n');

%% Plots -------------------------------------------------------------------
cols = [0 0.447 0.741; 0.85 0.325 0.098; 0.466 0.674 0.188];
for useNoise = [false true]
    figure('Name', sprintf('Closed loop, noise = %d', useNoise));
    sel = res([res.noise] == useNoise);
    panels = {1, 'Tank 1 height h_1 [cm]'; 2, 'Tank 2 height h_2 [cm]'; ...
              3, 'Tank 3 height h_3 [cm]'; 4, 'Tank 4 height h_4 [cm]'};
    for k = 1:4
        subplot(3, 2, k); hold on; grid on;
        for c = 1:numel(sel)
            plot(t/60, sel(c).h(panels{k,1}, :), 'Color', cols(c,:), 'LineWidth', 1.4);
        end
        if k <= 2, stairs(t/60, r(k,:), 'k--'); end
        ylabel(panels{k,2}); xlabel('Time [min]');
        if k == 1, legend({sel.name, 'reference'}, 'Location', 'southeast'); end
    end
    for j = 1:2
        subplot(3, 2, 4 + j); hold on; grid on;
        for c = 1:numel(sel)
            stairs(t/60, sel(c).u(j, :), 'Color', cols(c,:), 'LineWidth', 1.2);
        end
        ylabel(sprintf('Pump F_%d [cm^3/s]', j)); xlabel('Time [min]');
    end
end

%% ========================================================================
%  LOCAL FUNCTIONS (must be at the end of a MATLAB script)
% =========================================================================
function s = simulateClosedLoop(c, p, x0, t, r, u_bias, u_min, u_max, w, v)
    % Decentralised control: loop 1 = F1 -> h1, loop 2 = F2 -> h2.
    N = numel(t); Ts = t(2) - t(1);
    x = x0; h = zeros(4, N); u = zeros(2, N);
    st = struct('I', [0; 0], 'D', [0; 0], 'yprev', []);
    for k = 1:N
        h(:, k) = sys_sense(x, p);                       % true levels
        y = sys_out_z(x, p) + v(1:2, k);                 % measured h1, h2
        [u(:, k), st] = PIDControl(r(:, k), y, st, c, u_bias, u_min, u_max, Ts);
        if k < N                                         % ZOH over one sample
            [~, X] = ode15s(@(tt, xx) QuadrupleTankProcess(tt, xx, u(:,k) + w(:,k), p), ...
                            [t(k) t(k+1)], x);
            x = max(X(end, :).', 0);
        end
    end
    s = struct('name', c.name, 'h', h, 'u', u);
end

function [u, st] = PIDControl(r, y, st, c, u_bias, u_min, u_max, Ts)
    e = r - y;
    P = c.Kc * e;
    if c.Td > 0 && ~isempty(st.yprev)
        st.D = c.Tf/(c.Tf + Ts) * st.D - c.Kc*c.Td/(c.Tf + Ts) * (y - st.yprev);
    end
    st.yprev = y;
    v_unsat = u_bias + P + st.I + st.D;
    u = min(max(v_unsat, u_min), u_max);
    if isfinite(c.Ti)
        % integrate only if not saturated, or if the error drives u back inside
        ok = (v_unsat == u) | (sign(e) ~= sign(v_unsat - u));
        st.I = st.I + ok .* (c.Kc * Ts / c.Ti) .* e;
    end
end

function h = steadyStateLevels(F, p)
    a = p(1:4); gam = p(9:10); g = p(11);
    q3 = (1-gam(2))*F(2); q4 = (1-gam(1))*F(1);
    q  = [gam(1)*F(1) + q3; gam(2)*F(2) + q4; q3; q4];
    h  = q.^2 ./ (2*g*a.^2);
end

function [A, B, C] = linearise(h, p)
    a = p(1:4); At = p(5:8); gam = p(9:10); g = p(11); rho = p(12);
    T = At ./ a .* sqrt(2 * h / g);                     % tank time constants [s]
    A = diag(-1 ./ T); A(1,3) = 1/T(3); A(2,4) = 1/T(4);
    B = rho * [gam(1) 0; 0 gam(2); 0 1-gam(2); 1-gam(1) 0];
    C = diag(1 ./ (rho * At));
end

function z = tzero_2x2(A, B, C)
    n = size(A,1); m = size(B,2);
    M = [A B; C zeros(m)]; Nm = blkdiag(eye(n), zeros(m));
    z = eig(M, Nm); z = z(isfinite(z));
end

function xdot = QuadrupleTankProcess(~, x, u, p)
    m = x; F = u;
    a = p(1:4); A = p(5:8); gamma = p(9:10); g = p(11); rho = p(12);
    qin = [gamma(1)*F(1); gamma(2)*F(2); (1-gamma(2))*F(2); (1-gamma(1))*F(1)];
    h = max(m, 0) ./ (rho * A);          % clamp guards against tiny negative masses
    qout = a .* sqrt(2 * g * h);
    xdot = rho * [qin(1) + qout(3) - qout(1);
                  qin(2) + qout(4) - qout(2);
                  qin(3) - qout(3);
                  qin(4) - qout(4)];
end

function y = sys_sense(x, p)
    y = x ./ (p(12) * p(5:8));           % levels [cm]
end

function z = sys_out_z(x, p)
    y = sys_sense(x, p);
    z = y(1:2);                          % controlled variables h1, h2
end
