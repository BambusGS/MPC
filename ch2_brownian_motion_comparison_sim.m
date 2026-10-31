%% Disturbance Model Monte Carlo Comparison
% Compares three stochastic disturbance models for F3/F4:
% 1. Pure Brownian Motion (Random Walk)
% 2. Mean-Reverting Process (Ornstein-Uhlenbeck)
% 3. Square-Root Diffusion (keeps values non-negative by making diffusion dependent on the disturbance magnitude)

clear; clc; close all;

%% Configuration Parameters
n_runs = 8;          % Number of simulation runs to generate
T_total = 100;         % Total simulation time [seconds or minutes]
dt = 0.1;            % Time step size
t = 0:dt:T_total;
N = length(t);

% Model Parameters
F_bar = 50.0;          % Mean value
sigma = 2;          % Diffusion intensity 
a_rev = 0.1;          % Mean-reversion rate
F_init = 0.0;         % Initial starting value for all runs

% Preallocate storage matrices (rows = time steps, columns = runs)
F_brownian = zeros(N, n_runs);
F_mean_rev  = zeros(N, n_runs);
F_sqrt_diff = zeros(N, n_runs);

%% Sim (euler)
for j = 1:n_runs
    F_b = F_init;
    F_m = F_init;
    F_s = F_init;
    
    for i = 1:N       
        F_brownian(i, j) = F_b;
        F_mean_rev(i, j)  = F_m;
        F_sqrt_diff(i, j) = F_s;
        
        % Avoid stepping past the end of the time vector
        if i < N            
            dW = sqrt(dt) * randn(); 
            
            % Model 1: Pure Brownian Motion (Drift = 0, Diffusion = sigma)
            % dF = sigma * dw
            dF_b = sigma * dW;
            F_b = F_b + dF_b;
            
            % Model 2: Mean-Reverting
            % dF = a*(F_bar - F)dt + sigma * dw
            dF_m = a_rev * (F_bar - F_m) * dt + sigma * dW;
            F_m = F_m + dF_m;
            
            % Model 3: Square-Root State-Dependent Diffusion
            % dF = a*(F_bar - F)dt + sigma * sqrt(F) * dw          
            F_s_pos = max(F_s, 0); %sqrt non-negative fix
            dF_s = a_rev * (F_bar - F_s) * dt + sigma * sqrt(F_s_pos) * dW;
            F_s = F_s + dF_s;
        end
    end
end

%% Visualization

% Plot 1: Pure Brownian Motion
fig1 = figure('Name', 'Brownian Motion');
plot(t, F_brownian, 'LineWidth', 1.2);
grid on;
title({'Model 1: Brownian Motion', '$dF(t) = \sigma \, dw(t)$'}, 'Interpreter', 'latex', 'FontSize', 12);
xlabel('Time $t$ [s]', 'Interpreter', 'latex', 'FontSize', 11);
ylabel('Disturbance $F_3(t)$', 'Interpreter', 'latex', 'FontSize', 11);
xlim([0 T_total]);
legend('Location', 'best', 'Interpreter', 'latex');


% Plot 2: Mean-Reverting Process
fig2 = figure('Name', 'Mean Reverting Process');
plot(t, F_mean_rev, 'LineWidth', 1.2);
hold on;
plot(t, F_bar * ones(size(t)), 'k--', 'LineWidth', 2, 'DisplayName', 'Nominal Mean $\bar{F}_3$');
hold off;
grid on;
title({'Model 2: Mean-Reverting', '$dF(t) = a(\bar{F} - F(t))dt + \sigma \, dw(t)$'}, 'Interpreter', 'latex', 'FontSize', 12);
xlabel('Time $t$ [s]', 'Interpreter', 'latex', 'FontSize', 11);
ylabel('Disturbance $F_3(t)$', 'Interpreter', 'latex', 'FontSize', 11);
xlim([0 T_total]);
legend('Location', 'best', 'Interpreter', 'latex');


% Plot 3: Square-Root Diffusion Process
fig3 = figure('Name', 'Square Root Diffusion');
plot(t, F_sqrt_diff, 'LineWidth', 1.2);
hold on;
plot(t, F_bar * ones(size(t)), 'k--', 'LineWidth', 2, 'DisplayName', 'Nominal Mean $\bar{F}_3$');
hold off;
grid on;
title({'Model 3: Square-Root Diffusion (Poisson-type)', '$dF(t) = a(\bar{F} - F(t))dt + \sigma \sqrt{F(t)} \, dw(t)$'}, 'Interpreter', 'latex', 'FontSize', 12);
xlabel('Time $t$ [s]', 'Interpreter', 'latex', 'FontSize', 11);
ylabel('Disturbance $F_3(t)$', 'Interpreter', 'latex', 'FontSize', 11);
xlim([0 T_total]);
legend('Location', 'best', 'Interpreter', 'latex');

%% Export
disp('Plots generated and saved as PNG files:');
disp(' - figures/model1_brownian_motion.png');
exportgraphics(fig1, 'figures/model1_brownian_motion.png', 'Resolution', 600);
disp(' - figures/model2_mean_reverting.png');
exportgraphics(fig2, 'figures/model2_mean_reverting.png', 'Resolution', 600);
disp(' - figures/model3_sqrt_diffusion.png');
exportgraphics(fig3, 'figures/model3_sqrt_diffusion.png', 'Resolution', 600);