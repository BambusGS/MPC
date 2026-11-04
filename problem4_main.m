%% PROBLEM 4 - NONLINEAR SIMULATION
clear;
clc;
close all;

%% 1. ADD FUNCTION FOLDER TO MATLAB PATH

addpath("functions");

%% 2. CREATE SYSTEM PARAMETERS AND SIMULATION SETTINGS

p = createParameters();
sim = createSimSettings(p);

%% 3. DEFINE THE STEADY-STATE OPERATING POINT

% Pump flow rates [cm^3/s]
F1s = 300;
F2s = 300;

% Steady-state input vector
u_ss = [F1s; F2s];

% Deterministic

% For the deterministic baseline, no external disturbances
d_ss = [0; 0];

% Initial guess for the steady-state tank masses [g]
x_guess = 5000 * ones(4,1);

% Calculate steady state by solving dx/dt = 0
x_ss = fsolve(@(x) quadrupleTankProcess( ...
    0, x, u_ss, d_ss, p), x_guess);

% Calculate steady-state tank heights
y_ss = systemSense(x_ss, p);

% Stochastic

% For the Stochastic noise 
d_mean = p.FMean;

% Calculate steady state by solving dx/dt = 0
x_ss_noise = fsolve(@(x) quadrupleTankProcess( ...
    0, x, u_ss, d_mean, p), x_guess);

% Calculate steady-state tank heights
y_ss_noise = systemSense(x_ss_noise, p);




%% 4. INITIALIZE DETERMINISTIC SIMULATION

% Start every experiment from the calculated steady state
sim.initialState = x_ss;

% No external disturbances in the deterministic experiment
sim.initialDisturbance = d_ss;

% Constant baseline pump inputs
sim.input = repmat(u_ss, 1, sim.N);

%% 5. PART 4.1 - DETERMINISTIC STEP RESPONSES

step_sizes = [0.10, 0.25, 0.50];

% Store all simulation results
step_results = cell(2, length(step_sizes));

% Loop through both pumps
for pump = 1:2

    figure;
    tiledlayout(1,2);

    sgtitle(sprintf('Deterministic step responses - Pump %d', pump));

    % Loop through the three step sizes
    for i = 1:length(step_sizes)

        step_size = step_sizes(i);

        % Start from the steady-state
        u_step = u_ss;

        % Increase the selected pump
        u_step(pump) = u_ss(pump) * (1 + step_size);

        % Create a separate simulation setup
        sim_case = sim;

        % Apply the new constant pump inputs
        sim_case.input = repmat(u_step, 1, numel(sim.t));

        % Run the deterministic simulation
        result = simulateDeterministic(p, sim_case);

        % Store the result
        step_results{pump,i} = result;

        % Plot measured tank heights
        for tank = 1:2

            nexttile(tank);
            hold on;
            grid on;

            plot(result.t/60, result.y(:,tank), ...
                'LineWidth', 1.5, ...
                'DisplayName', ...
                sprintf('+%.0f%%', step_size*100));

            xlabel('Time [min]');
            ylabel(sprintf('Tank %d height [cm]', tank));
            title(sprintf('Tank %d', tank));
        end
    end

    % Add legends
    nexttile(1);
    legend('Location','best');

    nexttile(2);
    legend('Location','best');
end


%% 6. PART 4.2 - STOCHASTIC STEP RESPONSES

rng(68); % Reproducible random noise

% Three noise levels
noise_scaling = [0.5, 1.0, 2.0];
noise_name = ["Low", "Medium", "High"];

% make space for results
noisy_results = cell(2, 3, 3);

% Loop through noise levels
for noise_level = 1:3

    % Select current noise scaling
    scale = noise_scaling(noise_level);

    % Copy original parameters
    p_sto = p;

    % Scale process noise
    p_sto.piecewiseSigma = scale * p.piecewiseSigma;

    % Scale measurement noise
    p_sto.measurementCovariance = ...
        scale^2 * p.measurementCovariance;

    % Loop through both pumps
    for pump = 1:2

        figure;
        tiledlayout(1,2);

        sgtitle(sprintf('%s noise - Pump %d', ...
            noise_name(noise_level), pump));

        % Loop through step sizes
        for i = 1:length(step_sizes)

            step_size = step_sizes(i);

            % Define pump inputs
            u_step = u_ss;
            u_step(pump) = u_ss(pump)*(1 + step_size);

            % Configure simulation
            sim_case = sim;

            sim_case.initialState = x_ss_noise;
            sim_case.initialDisturbance = d_mean;

            sim_case.input = repmat(u_step, 1, sim.N);

            % Run stochastic simulation
            result = simulatePiecewise(p_sto, sim_case);

            % Store results
            noisy_results{pump,i,noise_level} = result;

            % Plot measured tank heights
            for tank = 1:2

                nexttile(tank);
                hold on;
                grid on;

                plot(result.t/60, result.y(:,tank), ...
                    'LineWidth', 1.2, ...
                    'DisplayName', ...
                    sprintf('+%.0f%%', step_size*100));

                xlabel('Time [min]');
                ylabel(sprintf('Tank %d height [cm]',tank));
                title(sprintf('Tank %d',tank));
            end
        end

        nexttile(1);
        legend('Location','best');

        nexttile(2);
        legend('Location','best');

    end
end

%% 7. PART 4.3 - NORMALIZED STEP RESPONSE

% space for results
norm_det = cell(2,3);
norm_noise = cell(2,3,3);

% Normalize deterministic responses
for pump = 1:2
    for i = 1:3

        result = step_results{pump,i};

        % Find steady state from the last 10%
        N = size(result.y,1);
        N_tail = max(1,round(0.1*N));

        y_initial = result.y(1,:);
        y_final = mean(result.y(end-N_tail+1:end,:),1);

        % Normalize
        norm_det{pump,i} = ...
            (result.y-y_initial)./(y_final-y_initial);

    end
end

% Normalize noisy responses
for pump = 1:2
    for i = 1:3
        for noise_level = 1:3

            result = noisy_results{pump,i,noise_level};

            % Find steady state from the last 10%
            N = size(result.y,1);
            N_tail = max(1,round(0.1*N));

            y_initial = result.y(1,:);
            y_final = mean(result.y(end-N_tail+1:end,:),1);

            % Normalize
            norm_noise{pump,i,noise_level} = ...
                (result.y-y_initial)./(y_final-y_initial);

        end
    end
end

        
% Plot determanistic

for pump = 1:2

    figure;
    tiledlayout(1,2);

    sgtitle(sprintf('Deterministic normalized responses - Pump %d',pump));

    for tank = 1:2

        nexttile;
        hold on;
        grid on;

        for i = 1:3

            result = step_results{pump,i};
            y_norm = norm_det{pump,i};

            plot(result.t/60,y_norm(:,tank), ...
                'LineWidth',1.5, ...
                'DisplayName',sprintf('+%.0f%%',step_sizes(i)*100));

        end

        xlabel('Time [min]');
        ylabel('Normalized tank height [-]');
        title(sprintf('Tank %d',tank));
        legend('Location','best');

    end
end


% Plot noisy

for pump = 1:2

    for noise_level = 1:3

        figure;
        tiledlayout(1,2);

        sgtitle(sprintf('Noisy normalized responses - Pump %d - %s', ...
            pump,noise_name(noise_level)));

        for tank = 1:2

            nexttile;
            hold on;
            grid on;

            for i = 1:3

                result = noisy_results{pump,i,noise_level};
                y_norm = norm_noise{pump,i,noise_level};

                plot(result.t/60,y_norm(:,tank), ...
                    'LineWidth',1.2, ...
                    'DisplayName',sprintf('+%.0f%%',step_sizes(i)*100));

            end

            xlabel('Time [min]');
            ylabel('Normalized tank height [-]');
            title(sprintf('Tank %d',tank));
            legend('Location','best');

        end
    end
end

%% 8. PART 4.4 - Tranfer-function

% Calculate steady-state gains
K = zeros(2,2,3);

for pump = 1:2
    for i = 1:3

        % run same code as in 4.3 for finding delta u and calculate gain
        result = step_results{pump,i};

     
        N = size(result.y,1);
        N_tail = max(1,round(0.1*N));

        y_initial = result.y(1,:);
        y_final = mean(result.y(end-N_tail+1:end,:),1);

        % Calculate input step
        delta_u = step_sizes(i)*u_ss(pump);

        % Calculate gains for both tanks
        K(:,pump,i) = ((y_final-y_initial)/delta_u).';

    end
end


% Average gains across the three step sizes
K_mean = mean(K,3);

disp('Identified steady-state gain matrix:');
disp(K_mean);


% Fit transfer functions

% Initial guesses for [T1, T2, beta, tau]
initial_guess = [100, 300, 10, 0];

% Lower and upper parameter limits
lower_bounds = [0.01, 0.01, 0, 0];
upper_bounds = [5000, 5000, 5000, 300];

% Store estimated parameters
estimated_params = zeros(2,2,4);

% Time vector in seconds
t = step_results{1,1}.t;

% Fit all four input-output relationships
for pump = 1:2
    for tank = 1:2

        % Average the three normalized step responses
        y_average = zeros(length(t),1);

        for i = 1:3
            y_average = y_average + norm_det{pump,i}(:,tank)/3;
        end

        % Fit transfer function
        params = lsqcurvefit( ...
            @transferStepResponse, ...
            initial_guess, ...
            t, ...
            y_average, ...
            lower_bounds, ...
            upper_bounds);

        % Store parameters
        estimated_params(tank,pump,:) = params;

    end
end

disp('Estimated transfer function parameters:');
disp(estimated_params);


function y = transferStepResponse(params,t)

    T1 = params(1);
    T2 = params(2);
    beta = params(3);
    tau = params(4);

    G = tf([beta 1],[T1*T2 T1+T2 1],'InputDelay',tau);

    y = squeeze(step(G,t));

end

% Compare fitted models with normalized responses

figure;
tiledlayout(2,2);

for pump = 1:2
    for tank = 1:2

        nexttile;
        hold on;
        grid on;

        % Get estimated parameters
        params = squeeze(estimated_params(tank,pump,:));

        % Calculate average measured response
        y_average = zeros(length(t),1);

        for i = 1:3
            y_average = y_average + norm_det{pump,i}(:,tank)/3;
        end

        % Calculate fitted response
        y_fit = transferStepResponse(params,t);

        % Plot both responses
        plot(t/60,y_average,'b','LineWidth',1.5);
        plot(t/60,y_fit,'r--','LineWidth',1.5);

        xlabel('Time [min]');
        ylabel('Normalized response');

        title(sprintf('Pump %d to Tank %d',pump,tank));

        legend('Simulated response','Transfer function');
    end
end

% Construct final transfer function matrix

G = cell(2,2);

for tank = 1:2
    for pump = 1:2

        T1 = estimated_params(tank,pump,1);
        T2 = estimated_params(tank,pump,2);
        beta = estimated_params(tank,pump,3);

        G{tank,pump} = K_mean(tank,pump) * ...
            tf([beta 1], [T1*T2, T1+T2, 1]);

    end
end

G_final = [G{1,1} G{1,2};
           G{2,1} G{2,2}];

G_final

% Evaluate transfer function fitting error

for pump = 1:2
    for tank = 1:2

        y_average = zeros(length(t),1);

        for i = 1:3
            y_average = y_average + norm_det{pump,i}(:,tank)/3;
        end

        params = squeeze(estimated_params(tank,pump,:));
        y_fit = transferStepResponse(params,t);

        %find root mean square error

        RMSE = sqrt(mean((y_average-y_fit).^2));

        fprintf('Pump %d -> Tank %d: RMSE = %.4f\n', ...
            pump, tank, RMSE);

    end
end
%they fit well

%% 9. PART 4.6 - Impulse response coefficients

% Discretization

Ts = sim.Ts;

G_discrete = c2d(G_final, Ts, 'zoh');

disp('Discrete-time transfer function:');
G_discrete

% Calculate Markov parameters

N = 50;                  % Number of coefficients
markov = zeros(2,2,N);

for tank = 1:2
    for pump = 1:2

        % Extract transfer function
        Gij = G_discrete(tank,pump);

        % Calculate impulse response
        t_impulse = (0:N-1) * Ts;
        h = impulse(Gij, t_impulse);

        % Store coefficients
        markov(tank,pump,:) = h;

    end
end

% Display Markov parameters

time = (0:N-1) * Ts;

figure;
tiledlayout(2,2);

for tank = 1:2
    for pump = 1:2

        nexttile;

        h = squeeze(markov(tank,pump,:));

        stem(time, h, 'filled');

        grid on;
        xlabel('Time [seconds]');
        ylabel('Markov parameter');
        title(sprintf('Pump %d -> Tank %d', pump, tank));

    end
end

sgtitle('Discrete-time impulse responses (Markov parameters)');