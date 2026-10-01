% Implement the four tank system in Matlab and simulate it without a controller
%   (open-loop simulation) and simulate it with PID controllers (closed-loop simulation).

clear;
clc;
close all;

%% PART1: non-linear model setup, and simulation ODE45
% --- PARAMETERS ---
p = struct();

p.a1 = 1.2272; % Pipe cross section for tank 1 [cm^2]
p.a2 = 1.2272; % Pipe cross section for tank 2 [cm^2]
p.a3 = 1.2272; % Pipe cross section for tank 3 [cm^2]
p.a4 = 1.2272; % Pipe cross section for tank 4 [cm^2]

p.A1 = 380.1327; % Tank cross sectional area for tank 1 [cm^2]
p.A2 = 380.1327; % Tank cross sectional area for tank 2 [cm^2]
p.A3 = 380.1327; % Tank cross sectional area for tank 3 [cm^2]
p.A4 = 380.1327; % Tank cross sectional area for tank 4 [cm^2]

p.gamma1 = 0.65; % Valve position for tank 1 [-]
p.gamma2 = 0.55; % Valve position for tank 2 [-]

p.g = 981; % gravity acceleration [cm/s^2]
p.rho = 1.00; % density of water [g/cm^3]

p_vector = [p.a1; p.a2; p.a3; p.a4; p.A1; p.A2; p.A3; p.A4; p.gamma1; p.gamma2; p.g; p.rho]

% --- SIM SETUP ---
t_0 = 0; % Start time [s]
t_f = 20*60; % End time [s]

T_s = 10; % Sampling time [s]        
t = [t_0:T_s:t_f]'; % Time vector [s]
N = length(t); % Number of time steps

m1_0 = 0; % Initial mass in tank 1 [g]
m2_0 = 0; % Initial mass in tank 2 [g]
m3_0 = 0; % Initial mass in tank 3 [g]
m4_0 = 0; % Initial mass in tank 4 [g]

F1 = 300; % Flow rate in pump 1 [cm^3/s]
F2 = 300; % Flow rate in pump 2 [cm^3/s]

x_0 = [m1_0; m2_0; m3_0; m4_0]; % Initial state vector
u = [repmat(F1, 1, N); repmat(F2, 1, N)]; % Input vector, repeated for N time steps
% -> CV z(controlled variables, subset of observed states Y) tank 1 and 2 heights, MV u(manipulated variables) are pump flow rates that are acting on CV
ref = [30; 30]; % Reference signal for tank heights 1 and 2 [cm]
Kp = 0.2; % Proportional gain for PI controller
Tinteg = 5000; % Integral time constant for PI controller [s]
integral_term = [0; 0]; % Initial integral term for PI controller

%% FUNCTION DEF
function xdot = QuadrupleTankProcess(t, x, u, p)
    % Four tank system model: dx/dt = f(t,x,u,p)
    %
    % This function implements a differential equation model
    % for the 4-tank system.

    % Unpack params
    m = x; % Mass of liquid in each tank [g]
    F = u; % Flow rates in pumps [cm3/s]
    a = p(1:4, 1); % Pipe cross sections [cm2]
    A = p(5:8, 1); % Tank cross sectional areas [cm2]
    gamma = p(9:10, 1); % Valve positions [-]
    g = p(11, 1); % Acceleration of gravity [cm/s2]
    rho = p(12, 1); % Density of water [g/cm3]

    % Inflows
    qin = zeros(4, 1);
    qin(1, 1) = gamma(1) * F(1); % Valve 1 to tank 1 [cm3/s]
    qin(2, 1) = gamma(2) * F(2); % Valve 2 to tank 2 [cm3/s]
    qin(3, 1) = (1 - gamma(2)) * F(2); % Valve 2 to tank 3 [cm3/s]
    qin(4, 1) = (1 - gamma(1)) * F(1); % Valve 1 to tank 4 [cm3/s]

    % Outflows
    h = m ./ (rho * A); % Liquid level in each tank [cm]
    qout = a .* sqrt(2 * g * h); % Outflow from each tank [cm3/s]

    % Differential equations for mass balances
    xdot = zeros(4, 1);
    xdot(1, 1) = rho * (qin(1, 1) + qout(3, 1) - qout(1, 1)); % Tank 1
    xdot(2, 1) = rho * (qin(2, 1) + qout(4, 1) - qout(2, 1)); % Tank 2
    xdot(3, 1) = rho * (qin(3, 1) - qout(3, 1)); % Tank 3
    xdot(4, 1) = rho * (qin(4, 1) - qout(4, 1)); % Tank 4
end

function [u, i] = PIControl(integral_term, ref, y, us, Kc, Ti, Ts)
    % PI controller implementation
    % Inputs:
    %   integral_term: integral term from previous step
    %   ref: reference signal (setpoint)
    %   y: measured output
    %   us: previous control signal
    %   Kc: proportional gain
    %   Ti: integral time constant
    %   Ts: sampling time
    %
    % Outputs:
    %   u: new control signal
    %   i: updated integral term

    e = ref - y; % Error signal
    u = us + Kc * e + integral_term; % Compute new control signal
    i = integral_term + (Kc * Ts / Ti) * e; % Update integral term
end


function y = sys_sense(x, p)
    % This function converts the state vector x into measured outputs y based on the system parameters p.
    % This system x is in grams, y is in cm, and p contains the parameters for the system.
    % Unpack params
    rho = p(12, 1); % Density of water [g/cm3]
    A = p(5:8, 1); % Tank cross sectional areas [cm2]
    % Convert state to measured outputs
    y = x ./ (rho * A); % Measured outputs (tank heights in cm)
end


function z = sys_out_z(x, p)
    % This function gets the controlled variables z from the observed states y.
    % For the four-tank system, the controlled variables are the heights of tanks 1 and 2.
    y = sys_sense(x, p);
    z = y(1:2); % Controlled variables (tank heights in cm)
end

function plotTankHeights(title, calcTime, calcHeights, measTime, measuredHeights)
    % Plot heights for all four tanks.
    % Heights must be N-by-4.
    % Measured data is optional.

    figure(1);
    clf(1);
    sgtitle(title);

    tankOrder = [3 4 1 2];
    colors = [
        0.0000 0.4470 0.7410
        0.8500 0.3250 0.0980
        0.9290 0.6940 0.1250
        0.4940 0.1840 0.5560
    ];

    hasMeasuredData = nargin >= 4 && ~isempty(measTime) && ~isempty(measuredHeights);

    for plotIndex = 1:4
        tankIndex = tankOrder(plotIndex);

        subplot(2, 2, plotIndex);
        hold on;

        plot(calcTime / 60, calcHeights(:, tankIndex), "-", "Color", colors(tankIndex, :), "LineWidth", 1.8);

        if hasMeasuredData
            plot(measTime / 60, measuredHeights(:, tankIndex), "--", "Color", colors(tankIndex, :), "LineWidth", 1.8);

            legend("Calculated", "Measured");
        end

        xlabel("Time [min]");
        ylabel(sprintf("Tank %d Height [cm]", tankIndex));
        grid on;
    end
end

function plotPumpFlows(title, time, pumpFlows)
    % Plot the two pump flow rates.
    % pumpFlows must be 2-by-N.

    figure(2);
    clf(2);
    sgtitle(title);

    colors = [
        0.0000 0.4470 0.7410
        0.8500 0.3250 0.0980
    ];

    for pumpIndex = 1:2
        subplot(2, 1, pumpIndex);

        plot(time / 60, pumpFlows(pumpIndex, :).', "Color", colors(pumpIndex, :), "LineWidth", 1.5);

        xlabel("Time [min]");
        ylabel(sprintf("Pump %d Flow Rate [cm^3/s]", pumpIndex));
        grid on;
    end
end


%% GAUSSIAN NOISE
seed = 100;
randn('state', seed); % Set random seed for reproducibility

% Process noise
Q = [20^2 0; 0 40^2]; % Covariance matrix for process noise
Lq = chol(Q, 'lower'); % Cholesky decomposition for generating correlated noise
w = Lq*randn(2, N); % Generate process noise for N time steps    

% Measurement noise
R = eye(4);
Lr = chol(R, 'lower'); % Cholesky decomposition for generating correlated noise
v = Lr*randn(4, N); % Generate measurement noise for N time


%% RUN SIM OPEN LOOP
% discrete time solver, with X representing the state vector over time
[T, X] = ode15s(@(t, x) QuadrupleTankProcess(t, x, u, p_vector), [t_0 t_f], x_0); % solver for stiff (sudden sharp changes) systems, try ode45 for non-stiff systems
% [outputs] = diff_eq_solver(@(t, x) function(t, x, parameters), [start_time end_time], initial_conditions)

% Helper variables to extract vector of tank heights and flow rates from the state vector
[nT,nX] = size(X);
a = p_vector(1:4,1)';
A = p_vector(5:8,1)';

% Compute measured variables (tank height in cm instead of X, which is mass in g)
H = zeros(nT,nX);
for i=1:nT
    H(i,:) = X(i,:)./(p.rho*A);
end

% Compute flows out from each tank
Qout = zeros(nT,nX);
for i=1:nT
    Qout(i,:) = a.*sqrt(2*p.g*H(i,:));
end

% Plot the tank heights over time in minutes
figure(1);
sgtitle("Tank Heights vs Time");

subplot(2,2,1);
plot(T/60, H(:,3), 'g', 'LineWidth', 1.5);
xlabel("Time [min]");
ylabel("Tank 3 Height [cm]");
grid on;

subplot(2,2,2);
plot(T/60, H(:,4), 'm', 'LineWidth', 1.5);
xlabel("Time [min]");
ylabel("Tank 4 Height [cm]");
grid on;

subplot(2,2,3);
plot(T/60, H(:,1), 'b', 'LineWidth', 1.5);
xlabel("Time [min]");
ylabel("Tank 1 Height [cm]");
grid on;

subplot(2,2,4);
plot(T/60, H(:,2), 'r', 'LineWidth', 1.5);
xlabel("Time [min]");
ylabel("Tank 2 Height [cm]");
grid on;

% Plot the pump flow rates over time in minutes
figure(2);
sgtitle("Pump Flow Rates vs Time");

subplot(2,1,1);
plot(T/60, u(1)*ones(nT,1), 'b', 'LineWidth', 1.5);
xlabel("Time [min]");
ylabel("Pump 1 Flow Rate [cm^3/s]");
grid on;

subplot(2,1,2);
plot(T/60, u(2)*ones(nT,1), 'r', 'LineWidth', 1.5);
xlabel("Time [min]");
ylabel("Pump 2 Flow Rate [cm^3/s]");
grid on;


%% RUN SIM CLOSED LOOP PI
% run this in discrete time steps

% define the number of states, inputs, outputs, and controlled variables and initialize the state, output, and controlled variable matrices
nx = 4; nu = 2; ny = 4; nz = 2;
x = zeros(nx,N);
y = zeros(ny,N);
u = zeros(nu,N);
z = zeros(nz,N);

X = zeros(0,nx);
T = zeros(0,1);

x(:,1) = x_0; % initial state
u(:,1) = [F1; F2]; % initial input
for i=1:N-1
    y(:, i) = sys_sense(x(:,i), p_vector) + v(:, i); % measured outputs + measurement noise
    z(:, i) = sys_out_z(x(:,i), p_vector) + v(1:2, i); % controlled variables + measurement noise

    % Controller
    [u(:, i+1), integral_term] = PIControl(integral_term, ref, z(:, i), u(:, i), Kp, Tinteg, T_s); % compute control action using PI controller

    % Model process
    [Ti, Xi] = ode15s(@(t, x) QuadrupleTankProcess(t, x, (u(:,i+1) + w(:,i)), p_vector), [t(i) t(i+1)], x(:,i)); % simulate system dynamics with process noise over one time step

    x(:,i+1) = Xi(end,:)'; % update state for next time step

    % Store the results for plotting
    X = [X; Xi]; % append the results of this time step to the overall state trajectory
    T = [T; Ti]; % append the time vector for this time step to the
end

i = N; % last time step, to complete the loop and get the final measurements
y(:, i) = sys_sense(x(:,i), p_vector) + v(:, i); % measured outputs + measurement noise
z(:, i) = sys_out_z(x(:,i), p_vector); % controlled variables


% Helper variables to extract vector of tank heights and flow rates from the state vector
[nT,nX] = size(X);
a = p_vector(1:4,1)';
A = p_vector(5:8,1)';

% Compute measured variables (tank height in cm instead of X, which is mass in g)
H = zeros(nT,nX);
for i=1:nT
    H(i,:) = X(i,:)./(p.rho*A);
end

H_meas = y';

% Compute flows out from each tank
Qout = zeros(nT,nX);
for i=1:nT
    Qout(i,:) = a.*sqrt(2*p.g*H(i,:));
end

% Plot the tank heights over time in minutes
figure(1);
clf(1);
sgtitle("Tank Heights vs Time");

% MATLAB default colors
c1 = [0.0000 0.4470 0.7410];   % Blue
c2 = [0.8500 0.3250 0.0980];   % Orange/red
c3 = [0.9290 0.6940 0.1250];   % Yellow
c4 = [0.4940 0.1840 0.5560];   % Purple

% Lighter colors for calculated values
c1_light = 0.5*c1 + 0.5*[1 1 1];
c2_light = 0.5*c2 + 0.5*[1 1 1];
c3_light = 0.5*c3 + 0.5*[1 1 1];
c4_light = 0.5*c4 + 0.5*[1 1 1];

subplot(2,2,1);
hold on;
plot(t/60, H_meas(:,3), '-', 'Color', c3, 'LineWidth', 1.8);
plot(T/60, H(:,3), '--', 'Color', c3_light, 'LineWidth', 1.8); 
xlabel("Time [min]");
ylabel("Tank 3 Height [cm]");
grid on;
legend("Calculated", "Measured");

subplot(2,2,2);
hold on;
plot(t/60, H_meas(:,4), '-', 'Color', c4, 'LineWidth', 1.8);
plot(T/60, H(:,4), '--', 'Color', c4_light, 'LineWidth', 1.8);
xlabel("Time [min]");
ylabel("Tank 4 Height [cm]");
grid on;
legend("Calculated", "Measured");

subplot(2,2,3);
hold on;
plot(t/60, H_meas(:,1), '-', 'Color', c1, 'LineWidth', 1.8);
plot(T/60, H(:,1), '--', 'Color', c1_light, 'LineWidth', 1.8);
xlabel("Time [min]");
ylabel("Tank 1 Height [cm]");
grid on;
legend("Calculated", "Measured");

subplot(2,2,4);
hold on;
plot(t/60, H_meas(:,2), '-', 'Color', c2, 'LineWidth', 1.8);
plot(T/60, H(:,2), '--', 'Color', c2_light, 'LineWidth', 1.8);
xlabel("Time [min]");
ylabel("Tank 2 Height [cm]");
grid on;
legend("Calculated", "Measured");

% Plot the pump flow rates over time in minutes
figure(2);
sgtitle("Pump Flow Rates vs Time");

subplot(2,1,1);
plot(t/60, u(1,:), 'b', 'LineWidth', 1.5);
xlabel("Time [min]");
ylabel("Pump 1 Flow Rate [cm^3/s]");
grid on;

subplot(2,1,2);
plot(t/60, u(2,:), 'r', 'LineWidth', 1.5);
xlabel("Time [min]");
ylabel("Pump 2 Flow Rate [cm^3/s]");
grid on;

%% SteadyState and Step Responses
% Consider the 4-tank system with the parameters given in p. 12. Let
% F1s = 250 cm3/s and F2s = 325 cm3/s.
% 1 Compute the steady-state of this system
% 2 Simulate the response for a 5, 10, and 25% step increase in F1,
% respectively.
% 3 Linearize the model at steady state. What is (A, B, C, D) for the
% continuous-time system?
% 4 Discretize the system using a sample time of Ts = 4 seconds. What is
% (A, B, C, D) for the discrete time system.
% 5 Simulate step responses for a 5, 10, and 25% step increase in F1,
% respectively using the linear model. Compare these responses to the
% responses of the nonlinear model.
% 6 Compute the continuous-time and discrete-time transfer functions for
% this 4-tank system. What are the gain, poles and zeros for these
% systems?

F1s = 250; % Steady-state flow rate for pump 1 [cm^3/s]
F2s = 325; % Steady-state flow rate for pump 2 [cm^3/s]

x_0 = 5000 * ones(4, 1); % Initial guess for steady-state mass in tanks [g]
u = [F1s; F2s]; % Steady-state input vector

x_ss = fsolve(@(x) QuadrupleTankProcess(0, x, u, p_vector), x_0); % Compute steady-state by solving dx/dt = 0
x_ss = real(x_ss); % Ensure steady-state is real-valued
y_ss = sys_sense(x_ss, p_vector); % Compute steady-state outputs (tank heights in cm)
z_ss = sys_out_z(x_ss, p_vector); % Compute steady-state controlled variables (tank heights in cm)

A = p_vector(5:8,1); % Tank cross sectional areas [cm^2]
H_ss = x_ss./(p.rho*A); % Steady-state heights in cm    
disp("Steady-state tank heights (cm):");
disp(H_ss);

%% sim and plot step responses
step_sizes = [0.05, 0.10, 0.25]; % Step sizes for F1 increase
for step_size = step_sizes
    F1_step = F1s * (1 + step_size); % New flow rate for pump 1 after step increase
    u_step = [F1_step; F2s]; % Input vector after step increase

    % Simulate the system response to the step input
    [T_step, X_step] = ode15s(@(t, x) QuadrupleTankProcess(t, x, u_step, p_vector), [t_0 t_f], x_ss);

    % Compute tank heights from mass
    nx = 4; % Number of tanks
    A = p_vector(5:8,1)'; % Tank cross sectional areas [cm^2]
    H_step = zeros(length(T_step), nx);
    for i = 1:length(T_step)
        H_step(i,:) = X_step(i,:)./(p.rho*A_lAin);
    end

    % plot the tank heights 3, 4, 1, 2 in cm over minutes, all sizes in one figure
    figure(3);
    clf(3);
    sgtitle("Tank Heights for F1 Step Increases");

    tankOrder = [3 4 1 2];

    for plotIndex = 1:4
        subplot(2, 2, plotIndex);
        hold on;
        grid on;

        tankIndex = tankOrder(plotIndex);

        for step_size = step_sizes
            F1_step = F1s * (1 + step_size);
            u_step = [F1_step; F2s];

            [T_step, X_step] = ode15s( ...
                @(t, x) QuadrupleTankProcess(t, x, u_step, p_vector), ...
                [t_0 t_f], ...
                x_ss);

            A = p_vector(5:8, 1)';
            H_step = X_step ./ (p.rho * A);

            plot( ...
                T_step / 60, ...
                H_step(:, tankIndex), ...
                "LineWidth", 1.5, ...
                "DisplayName", sprintf("F1 + %.0f%%", 100 * step_size));
        end

        xlabel("Time [min]");
        ylabel(sprintf("Tank %d Height [cm]", tankIndex));
        title(sprintf("Tank %d", tankIndex));
        legend("Location", "best");
    end
end

%% linearize the model around x_ss to get A, B, C, D matrices for continuous-time system
h_ss = y_ss; % Steady-state heights in cm

Area_tank = p_vector(5:8,1); % Tank cross sectional areas [cm^2]
area_pipe = p_vector(1:4,1); % Pipe cross sections [cm^2]
grav_accel = p_vector(11,1); % Gravity acceleration [cm/s^2]
rho = p_vector(12,1); % Density of water [g/cm^3]
gamma = p_vector(9:10,1); % Valve positions [-]

% obtained as Jacobian (manual derivatives for linearization)
T = Area_tank ./ area_pipe .* sqrt(2 * h_ss / grav_accel); % Time constants for each tank [s]; computed as A_tank / a_pipe * sqrt(2 * g * h_ss) [area ratio * gravity liquid level pull]

A_lin = diag(-1 ./ T);
A_lin(1, 3) = 1 / T(3); % Tank 1 affected by Tank 3 outflow
A_lin(2, 4) = 1 / T(4); % Tank 2 affected by Tank 4 outflow
B_lin = [rho*gamma(1), 0; 0, rho*gamma(2); 0, rho*(1-gamma(2)); rho*(1-gamma(1)), 0]; % Input matrix for pumps as derivative of dx/dt with respect to u (which is the pumps F1 and F2)
C = diag(1 ./ (rho .* Area_tank));
C_z = C(1:2, :); % Controlled variables are tank heights 1 and 2
D = zeros(4, 2); % No direct feedthrough from inputs to outputs
A_lin
B_lin
C


%% Discretize the system using a sample time of Ts = 4 seconds
Ts = 4; % Sample time [s]

% discretization happens through ZOH on the A_d and B_d matrices (Cd = C, Dd = D)
function [Abar,Bbar]=c2dzoh(A,B,Ts)
    [nx, nu] = size(B);
    M = [A, B; zeros(nu, nx + nu)];
    phi = expm(M * Ts);
    Abar = phi(1:nx, 1:nx);
    Bbar = phi(1:nx, nx + 1:end);
end

[A_d, B_d] = c2dzoh(A_lin, B_lin, Ts);
C_d = C;
D_d = D;

A_d
B_d
C_d 

%% plot the step responses for the discrete-time system for the same step sizes as before (dashed line of same color over the non-linear model response)
% has to be simulated in steps as delta_mass (heigh) that is then put on top of the x_ss
figure(3);
hold on;

tankOrder = [3 4 1 2];
stepColors = lines(length(step_sizes));

t_discrete = (t_0:Ts:t_f)';
N_discrete = length(t_discrete);

for stepIndex = 1:length(step_sizes)
    step_size = step_sizes(stepIndex);

    F1_step = F1s * (1 + step_size);
    delta_u = [F1_step; F2s] - [F1s; F2s];

    % Simulate deviation from steady state
    delta_x = zeros(4, N_discrete);
    for k = 1:N_discrete-1
        delta_x(:, k+1) = A_d * delta_x(:, k) + B_d * delta_u;
    end

    % Convert deviations to absolute tank heights
    H_discrete = (y_ss + C_d * delta_x).';

    % Overlay discrete-time response
    for plotIndex = 1:4
        subplot(2, 2, plotIndex);
        hold on;

        tankIndex = tankOrder(plotIndex);

        plot(t_discrete / 60, ...
            H_discrete(:, tankIndex), ...
            "--", ...
            "Color", stepColors(stepIndex, :), ...
            "LineWidth", 1.5, ...
            "DisplayName", sprintf("F1 + %.0f%% discrete", 100 * step_size));
    end
end

for plotIndex = 1:4
    subplot(2, 2, plotIndex);
    legend("Location", "best");
end

%% Compute the continuous-time and discrete-time transfer functions for this 4-tank system. What are the gain, poles and zeros for these systems
% basically Linear Ctrl II here
sys_ct = ss(A_lin, B_lin, C, D); 
G_ct = tf(sys_ct)

% poles(eigval) and zeroes(generalized eigval)
poles_ct = eig(A_lin) % the time constants T computed earlier basically
zeros_ct = tzero(sys_ct)

% DC gain
gain_ct = dcgain(sys_ct)


sys_d = ss(A_d, B_d, C_d, D_d, Ts);
G_dt = tf(sys_d)

poles_dt = eig(A_d) % the time constants basically
zeros_dt = tzero(sys_d)

gain_dt = dcgain(sys_d) % [!NOTE] is the same as the continuous-time system gain, as expected for a proper discretization

gain_ct - gain_dt % should be close to zero, as expected for a proper discretization 