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

% Constant external disturbances at their predefined mean values
d_ss = p.FMean;

% Initial guess for the steady-state tank masses [g]
x_guess = 5000 * ones(4,1);

% Calculate steady state by solving dx/dt = 0
x_ss = fsolve(@(x) quadrupleTankProcess( ...
    0, x, u_ss, d_ss, p), x_guess);

% Calculate steady-state tank heights
y_ss = systemSense(x_ss, p);

% Stochastic

x_ss_noise = x_ss;
y_ss_noise = y_ss;





%% 4. INITIALIZE DETERMINISTIC SIMULATION

% Start every experiment from the calculated steady state
sim.initialState = x_ss;

% Keep external disturbances constant at their mean values
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

    sgtitle(sprintf('Deterministic step responses - Pump F_%d', pump));

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
            ylabel(sprintf('Tank level h_%d [cm]', tank));
            title(sprintf('Pump F_%d to Tank h_%d', pump, tank));
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

        sgtitle(sprintf('%s-noise step responses - Pump F_%d', ...
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
            sim_case.initialDisturbance = d_ss;

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
                ylabel(sprintf('Tank level h_%d [cm]',tank));
                title(sprintf('Pump F_%d to Tank h_%d',pump,tank));
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

        y_initial = y_ss(:)';
        delta_u = step_sizes(i)*u_ss(pump);

        % Normalize output deviations by the corresponding pump input change
        norm_det{pump,i} = ...
            (result.y-y_initial) ./ delta_u;
    end
end

% Normalize noisy responses
for pump = 1:2
    for i = 1:3
        for noise_level = 1:3

            result = noisy_results{pump,i,noise_level};

            y_initial = y_ss_noise(:)';
            delta_u = step_sizes(i)*u_ss(pump);

            % Normalize
            norm_noise{pump,i,noise_level} = ...
                (result.y-y_initial) ./ delta_u;

        end
    end
end

        
% Plot deterministic

for pump = 1:2

    figure;
    tiledlayout(1,2);

    sgtitle(sprintf('Deterministic normalized step responses - Pump F_%d',pump));

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
        ylabel('Normalized response [cm/(cm^3/s)]');
        title(sprintf('Pump F_%d to Tank h_%d',pump,tank));
        legend('Location','best');

    end
end


% Plot noisy

for pump = 1:2

    for noise_level = 1:3

        figure;
        tiledlayout(1,2);

        sgtitle(sprintf('%s-noise normalized step responses - Pump F_%d', ...
            noise_name(noise_level),pump));

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
            ylabel('Normalized response [cm/(cm^3/s)]');
            title(sprintf('Pump F_%d to Tank h_%d',pump,tank));
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

        y_initial = y_ss(:)';
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


% Fit transfer functions using the structure from the lecture

% Parameter order:
% [K11, K12, K21, K22, T1, T2, T3, T4]

estimated_params_shared = zeros(8,3);

% Time vector in seconds
t = step_results{1,1}.t;

% Initial parameter estimates
initial_guess = [0.14, 0.068, 0.12, 0.22, 85, 115, 40, 55];

% Lower and upper parameter bounds
lower_bounds = [0, 0, 0, 0, 0.01, 0.01, 0.01, 0.01];
upper_bounds = [10, 10, 10, 10, 5000, 5000, 5000, 5000];

% Optimization settings
options = optimoptions('lsqcurvefit', ...
    'FunctionTolerance',1e-12, ...
    'StepTolerance',1e-12, ...
    'OptimalityTolerance',1e-10, ...
    'MaxIterations',2000, ...
    'MaxFunctionEvaluations',10000);

% Fit all four input-output responses simultaneously
for i = 1:3

    % Collect the four normalized responses:
    % Column 1: Pump 1 -> Tank 1 (G11)
    % Column 2: Pump 1 -> Tank 2 (G21)
    % Column 3: Pump 2 -> Tank 1 (G12)
    % Column 4: Pump 2 -> Tank 2 (G22)

    y_data = [ ...
        norm_det{1,i}(:,1), ...
        norm_det{1,i}(:,2), ...
        norm_det{2,i}(:,1), ...
        norm_det{2,i}(:,2)];

    % Estimate all eight parameters together
    params = lsqcurvefit( ...
        @sharedTransferStepResponse, ...
        initial_guess, ...
        t, ...
        y_data, ...
        lower_bounds, ...
        upper_bounds, ...
        options);

    % Store parameters for this step size
    estimated_params_shared(:,i) = params(:);

end

disp('Estimated parameters with shared time constants:');
disp(estimated_params_shared);
for i = 1:3
    params = estimated_params_shared(:,i);

    fprintf('\nEstimated parameters for %.0f%% step:\n', step_sizes(i)*100);
    fprintf('K11 = %.6g, K12 = %.6g, K21 = %.6g, K22 = %.6g\n', params(1:4));
    fprintf('T1  = %.6g, T2  = %.6g, T3  = %.6g, T4  = %.6g\n', params(5:8));
end


% PLOT THE SHARED MODEL AGAINST THE DETERMINISTIC RESPONSES

for i = 1:3

    params = estimated_params_shared(:,i);

    % Calculate all four fitted responses
    y_fit = sharedTransferStepResponse(params,t);

    figure;
    tiledlayout(2,2);

    sgtitle(sprintf('Deterministic transfer-function model - %.0f%% step', ...
        step_sizes(i)*100));

    % Pump 1 -> Tank 1
    nexttile;
    plot(t/60,norm_det{1,i}(:,1),'LineWidth',1.5);
    hold on;
    plot(t/60,y_fit(:,1),'--','LineWidth',1.5);
    grid on;
    xlabel('Time [min]');
    ylabel('Normalized response [cm/(cm^3/s)]');
    title('Pump F_1 to Tank h_1');
    legend('Simulated','Shared model','Location','best');

    % Pump 1 -> Tank 2
    nexttile;
    plot(t/60,norm_det{1,i}(:,2),'LineWidth',1.5);
    hold on;
    plot(t/60,y_fit(:,2),'--','LineWidth',1.5);
    grid on;
    xlabel('Time [min]');
    ylabel('Normalized response [cm/(cm^3/s)]');
    title('Pump F_1 to Tank h_2');
    legend('Simulated','Shared model','Location','best');

    % Pump 2 -> Tank 1
    nexttile;
    plot(t/60,norm_det{2,i}(:,1),'LineWidth',1.5);
    hold on;
    plot(t/60,y_fit(:,3),'--','LineWidth',1.5);
    grid on;
    xlabel('Time [min]');
    ylabel('Normalized response [cm/(cm^3/s)]');
    title('Pump F_2 to Tank h_1');
    legend('Simulated','Shared model','Location','best');

    % Pump 2 -> Tank 2
    nexttile;
    plot(t/60,norm_det{2,i}(:,2),'LineWidth',1.5);
    hold on;
    plot(t/60,y_fit(:,4),'--','LineWidth',1.5);
    grid on;
    xlabel('Time [min]');
    ylabel('Normalized response [cm/(cm^3/s)]');
    title('Pump F_2 to Tank h_2');
    legend('Simulated','Shared model','Location','best');

end


% Fit a new model using medium-noise stochastic data

pump = 1;
noise_level = 2;  % Medium noise

% Parameter order:
% [K11_10, K21_10, K11_25, K21_25, K11_50, K21_50, T1, T2, T4]

estimated_params_stochastic = zeros(9,1);

% Collect all six stochastic responses:
% Columns 1-2: 10% step, Tank 1 and Tank 2
% Columns 3-4: 25% step, Tank 1 and Tank 2
% Columns 5-6: 50% step, Tank 1 and Tank 2

y_stochastic_data = [];

for i = 1:3

    y_stochastic_data = [y_stochastic_data, ...
        norm_noise{pump,i,noise_level}(:,1), ...
        norm_noise{pump,i,noise_level}(:,2)];

end

% Use deterministic estimates as initial guesses for the gains
for i = 1:3
    estimated_params_stochastic(2*i-1) = estimated_params_shared(1,i);
    estimated_params_stochastic(2*i) = estimated_params_shared(3,i);
end

% Use average deterministic time constants as initial guesses
estimated_params_stochastic(7) = mean(estimated_params_shared(5,:));
estimated_params_stochastic(8) = mean(estimated_params_shared(6,:));
estimated_params_stochastic(9) = mean(estimated_params_shared(8,:));

% Parameter bounds
lower_bounds = [0, 0, 0, 0, 0, 0, 0.01, 0.01, 0.01];
upper_bounds = [10, 10, 10, 10, 10, 10, 5000, 5000, 5000];

% Optimization settings
options = optimoptions('lsqcurvefit', ...
    'FunctionTolerance',1e-12, ...
    'StepTolerance',1e-12, ...
    'OptimalityTolerance',1e-10, ...
    'MaxIterations',2000, ...
    'MaxFunctionEvaluations',10000);

% Fit all six responses simultaneously using shared time constants
estimated_params_stochastic = lsqcurvefit( ...
    @stochasticSharedStepResponse, ...
    estimated_params_stochastic, ...
    t, ...
    y_stochastic_data, ...
    lower_bounds, ...
    upper_bounds, ...
    options);

disp('Stochastic model parameters with shared time constants:');
disp(estimated_params_stochastic);


% PLOT STOCHASTIC DATA AND SHARED FITTED MODELS

figure;
tiledlayout(1,2);

sgtitle('Medium-noise stochastic responses and shared models - Pump F_1');

for tank = 1:2

    nexttile;
    hold on;
    grid on;

    for i = 1:3

        % Extract stochastic measurements
        result = noisy_results{pump,i,noise_level};
        y_stochastic = norm_noise{pump,i,noise_level}(:,tank);

        % Calculate fitted stochastic model response
        y_fit_stochastic = stochasticSharedStepResponse( ...
            estimated_params_stochastic,t);

        % Select the corresponding fitted response
        y_fit = y_fit_stochastic(:,2*i-2+tank);

        % Plot stochastic measurements
        plot(result.t/60,y_stochastic, ...
            'LineWidth',1.0, ...
            'DisplayName',sprintf('Stochastic data +%.0f%%', ...
            step_sizes(i)*100));

        % Plot fitted model
        plot(t/60,y_fit,'--', ...
            'LineWidth',1.5, ...
            'DisplayName',sprintf('Shared model +%.0f%%', ...
            step_sizes(i)*100));

    end

    xlabel('Time [min]');
    ylabel('Normalized response [cm/(cm^3/s)]');
    title(sprintf('Pump F_%d to Tank h_%d',pump,tank));
    legend('Location','best');

end

%% CONSTRUCT TRANSFER-FUNCTION MATRICES USING SHARED TIME CONSTANTS

G_final = cell(1,3);

for i = 1:3

    % Extract the eight estimated parameters
    params = estimated_params_shared(:,i);

    K11 = params(1);
    K12 = params(2);
    K21 = params(3);
    K22 = params(4);

    T1 = params(5);
    T2 = params(6);
    T3 = params(7);
    T4 = params(8);

    % Construct the transfer-function matrix
    G11 = tf(K11,[T1 1]);

    G12 = tf(K12,conv([T1 1],[T3 1]));

    G21 = tf(K21,conv([T2 1],[T4 1]));

    G22 = tf(K22,[T2 1]);

    G_final{i} = [G11 G12;
                  G21 G22];

end

disp('Transfer function matrix for 10% steps:');
G_final{1}

disp('Transfer function matrix for 25% steps:');
G_final{2}

disp('Transfer function matrix for 50% steps:');
G_final{3}


% CALCULATE RMSE FOR THE SHARED TRANSFER-FUNCTION MODELS

for i = 1:3

    % Extract parameters for this step size
    params = estimated_params_shared(:,i);

    % Calculate all four fitted responses
    y_fit = sharedTransferStepResponse(params,t);

    % Calculate RMSE for each input-output relationship

    % Pump 1 -> Tank 1
    y_data = norm_det{1,i}(:,1);
    RMSE = sqrt(mean((y_data-y_fit(:,1)).^2));
    fprintf('Pump 1 -> Tank 1, %.0f%% step: RMSE = %.6f\n', ...
        step_sizes(i)*100,RMSE);

    % Pump 1 -> Tank 2
    y_data = norm_det{1,i}(:,2);
    RMSE = sqrt(mean((y_data-y_fit(:,2)).^2));
    fprintf('Pump 1 -> Tank 2, %.0f%% step: RMSE = %.6f\n', ...
        step_sizes(i)*100,RMSE);

    % Pump 2 -> Tank 1
    y_data = norm_det{2,i}(:,1);
    RMSE = sqrt(mean((y_data-y_fit(:,3)).^2));
    fprintf('Pump 2 -> Tank 1, %.0f%% step: RMSE = %.6f\n', ...
        step_sizes(i)*100,RMSE);

    % Pump 2 -> Tank 2
    y_data = norm_det{2,i}(:,2);
    RMSE = sqrt(mean((y_data-y_fit(:,4)).^2));
    fprintf('Pump 2 -> Tank 2, %.0f%% step: RMSE = %.6f\n', ...
        step_sizes(i)*100,RMSE);

end

%% 9. PART 4.6 - Impulse response coefficients (Markov parameters)

% Discretization
Ts = sim.Ts;

% Discretize each transfer function matrix
G_discrete = cell(1,3);

for i = 1:3

    G_discrete{i} = c2d(G_final{i}, Ts, 'zoh');

    fprintf('Discrete-time transfer function for %.0f%% step:\n', ...
        step_sizes(i)*100);

    G_discrete{i}

end


%% Calculate Markov parameters

N = 90;  % Number of coefficients

% Store Markov parameters for all three step sizes
markov = zeros(2,2,N,3);

for i = 1:3

    for tank = 1:2
        for pump = 1:2

            Gij = G_discrete{i}(tank,pump);

            sys = ss(Gij);

            [A,B,C,D] = ssdata(sys);

            % First Markov parameter is the direct-feedthrough matrix
            markov(tank,pump,1,i) = D;

            % H1, H2, H3, ...
            for k = 2:N
                markov(tank,pump,k,i) = C*A^(k-2)*B;
            end

        end
    end
end

% save it for later use in problem 6
save("markov.mat","markov")


%% Plot Markov parameters for each step size

time = (0:N-1)*Ts;

for i = 1:3

    figure;
    tiledlayout(2,2);

    sgtitle(sprintf('Markov parameters - %.0f%% input step', ...
        step_sizes(i)*100));

    for tank = 1:2
        for pump = 1:2

            nexttile;

            h = squeeze(markov(tank,pump,:,i));

            stem(time,h,'filled');

            grid on;
            xlabel('Time [seconds]');
            ylabel('Markov parameter [cm/(cm^3/s)]');

            title(sprintf('Pump F_%d to Tank h_%d',pump,tank));

        end
    end
end




function y = stochasticSharedStepResponse(params,t)

    % Extract shared time constants
    T1 = params(7);
    T2 = params(8);
    T4 = params(9);

    % Initialize output matrix
    % Columns: G11_10, G21_10, G11_25, G21_25, G11_50, G21_50
    y = zeros(length(t),6);

    % Calculate responses for all three step sizes
    for i = 1:3

        % Extract gains for this step size
        K11 = params(2*i-1);
        K21 = params(2*i);

        % Pump 1 -> Tank 1 (first order)
        G11 = tf(K11,[T1 1]);

        % Pump 1 -> Tank 2 (second order)
        G21 = tf(K21,conv([T2 1],[T4 1]));

        % Calculate step responses
        y(:,2*i-1) = squeeze(step(G11,t));
        y(:,2*i) = squeeze(step(G21,t));

    end

end


function y = sharedTransferStepResponse(params,t)

    % Extract gains
    K11 = params(1);
    K12 = params(2);
    K21 = params(3);
    K22 = params(4);

    % Extract shared time constants
    T1 = params(5);
    T2 = params(6);
    T3 = params(7);
    T4 = params(8);

    % Construct the four transfer functions
    G11 = tf(K11,[T1 1]);

    G12 = tf(K12,conv([T1 1],[T3 1]));

    G21 = tf(K21,conv([T2 1],[T4 1]));

    G22 = tf(K22,[T2 1]);

    % Calculate step responses
    y = [ ...
        squeeze(step(G11,t)), ...
        squeeze(step(G21,t)), ...
        squeeze(step(G12,t)), ...
        squeeze(step(G22,t))];

end