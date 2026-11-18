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

        
% Plot determanistic

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

% Store estimated parameters:
% Diagonal:     [K, T]
% Off-diagonal: [K, Tmain, Textra]
estimated_params = zeros(2,2,3,3);

% Time vector in seconds
t = step_results{1,1}.t;

% Fit all four input-output relationships
for pump = 1:2
    for tank = 1:2
        for i = 1:3

            y_fit_data = norm_det{pump,i}(:,tank);

            options = optimoptions('lsqcurvefit', ...
                'FunctionTolerance',1e-12, ...
                'StepTolerance',1e-12, ...
                'OptimalityTolerance',1e-10, ...
                'MaxIterations',2000, ...
                'MaxFunctionEvaluations',10000);

            if tank == pump
                % Diagonal elements: first order [K, T]
                initial_guess = [0.1, 100];
                lower_bounds = [0, 0.01];
                upper_bounds = [10, 5000];
            else
                % Off-diagonal elements: second order [K, Tmain, Textra]
                initial_guess = [0.1, 100, 30];
                lower_bounds = [0, 0.01, 0.01];
                upper_bounds = [10, 5000, 5000];
            end

            % Fit using the correct model order
            params = lsqcurvefit( ...
                @(params,t) transferStepResponse(params,t,tank,pump), ...
                initial_guess, ...
                t, ...
                y_fit_data, ...
                lower_bounds, ...
                upper_bounds, ...
                options);
            % Store parameters for the current input-output pair and step size
            estimated_params(tank,pump,1:length(params),i) = params;

        end
    end
end

disp('Estimated transfer function parameters:');
disp(estimated_params);



% Compare fitted models for each step size

for pump = 1:2
    for tank = 1:2

        figure;
        hold on;
        grid on;

        for i = 1:3

            % Extract parameters for this step size
            params = squeeze(estimated_params(tank,pump,:,i));

            % Get corresponding simulated response
            result = step_results{pump,i};
            y_data = norm_det{pump,i}(:,tank);

            % Calculate fitted response
            y_fit = transferStepResponse(params,t,tank,pump);

            % Plot both responses
            plot(t/60,y_data,'LineWidth',1.5, ...
                'DisplayName',sprintf('Simulated +%.0f%%',step_sizes(i)*100));

            plot(t/60,y_fit,'--','LineWidth',1.5, ...
                'DisplayName',sprintf('Fitted model +%.0f%%',step_sizes(i)*100));

        end

        xlabel('Time [min]');
        ylabel('Normalized response [cm/(cm^3/s)]');
        title(sprintf('Pump F_%d to Tank h_%d',pump,tank));
        legend('Location','best');

    end
end



%% FIT A NEW MODEL USING MEDIUM-NOISE STOCHASTIC DATA

pump = 1;
noise_level = 2;  % Medium noise

% Store parameters for the stochastic models
% Diagonal: [K, T]
% Off-diagonal: [K, Tmain, Textra]
estimated_params_stochastic = zeros(2,2,3,3);

% Fit both tank responses for all three step sizes
for tank = 1:2
    for i = 1:3

        % Extract medium-noise stochastic response
        y_stochastic = norm_noise{pump,i,noise_level}(:,tank);

        % Select model structure and parameter limits
        if tank == pump
            initial_guess = [0.1, 100];
            lower_bounds = [0, 0.01];
            upper_bounds = [10, 5000];
        else
            initial_guess = [0.1, 100, 30];
            lower_bounds = [0, 0.01, 0.01];
            upper_bounds = [10, 5000, 5000];
        end

        % Configure least-squares optimization
        options = optimoptions('lsqcurvefit', ...
            'FunctionTolerance',1e-12, ...
            'StepTolerance',1e-12, ...
            'OptimalityTolerance',1e-10, ...
            'MaxIterations',2000, ...
            'MaxFunctionEvaluations',10000);

        % Fit model directly to stochastic data
        params_stochastic = lsqcurvefit( ...
            @(params,t) transferStepResponse(params,t,tank,pump), ...
            initial_guess, ...
            t, ...
            y_stochastic, ...
            lower_bounds, ...
            upper_bounds, ...
            options);

        % Store estimated parameters
        estimated_params_stochastic(tank,pump, ...
            1:length(params_stochastic),i) = params_stochastic;

    end
end

disp('Transfer-function parameters identified from medium-noise data:');
disp(estimated_params_stochastic);


%% PLOT STOCHASTIC DATA AND THEIR NEWLY FITTED MODELS

figure;
tiledlayout(1,2);

sgtitle('Medium-noise stochastic responses and identified models - Pump F_1');

for tank = 1:2

    nexttile;
    hold on;
    grid on;

    for i = 1:3

        % Extract stochastic measurements
        result = noisy_results{pump,i,noise_level};
        y_stochastic = norm_noise{pump,i,noise_level}(:,tank);

        % Extract parameters fitted to stochastic data
        params_stochastic = squeeze( ...
            estimated_params_stochastic(tank,pump,:,i));

        % Calculate newly fitted model response
        y_fit_stochastic = transferStepResponse( ...
            params_stochastic,t,tank,pump);

        % Plot stochastic measurements
        plot(result.t/60,y_stochastic, ...
            'LineWidth',1.0, ...
            'DisplayName',sprintf('Stochastic data +%.0f%%', ...
            step_sizes(i)*100));

        % Plot newly fitted model
        plot(t/60,y_fit_stochastic,'--', ...
            'LineWidth',1.5, ...
            'DisplayName',sprintf('Stochastic fitted model +%.0f%%', ...
            step_sizes(i)*100));

    end

    xlabel('Time [min]');
    ylabel('Normalized response [cm/(cm^3/s)]');
    title(sprintf('Pump F_%d to Tank h_%d',pump,tank));
    legend('Location','best');

end

% Construct transfer-function matrices for each step size

G_final = cell(1,3);

for i = 1:3

    % Create a new 2x2 matrix for this step size
    Gcell = cell(2,2);

    for tank = 1:2
        for pump = 1:2

            % Extract parameters for this specific input-output pair
            params = reshape(estimated_params(tank,pump,:,i),[],1);

            Kij = params(1);
            Tmain = params(2);

            if tank == pump
                % Diagonal elements: first order
                Gcell{tank,pump} = tf(Kij,[Tmain 1]);

            else
                % Off-diagonal elements: second order
                Textra = params(3);

                denominator = conv([Tmain 1],[Textra 1]);

                Gcell{tank,pump} = tf(Kij,denominator);
            end

        end
    end

    % Assemble the 2x2 transfer-function matrix
    G_final{i} = [Gcell{1,1}, Gcell{1,2};
                  Gcell{2,1}, Gcell{2,2}];

end

disp('Transfer function matrix for 10% steps:');
G_final{1}

disp('Transfer function matrix for 25% steps:');
G_final{2}

disp('Transfer function matrix for 50% steps:');
G_final{3}

% Calculate RMSE for each step size

for pump = 1:2
    for tank = 1:2
        for i = 1:3

            % Extract fitted parameters
            params = squeeze(estimated_params(tank,pump,:,i));

            % Get simulated response
            y_data = norm_det{pump,i}(:,tank);

            % Calculate fitted response
            y_fit = transferStepResponse(params,t,tank,pump);

            % Calculate RMSE
            RMSE = sqrt(mean((y_data-y_fit).^2));

            fprintf('Pump %d -> Tank %d, %.0f%% step: RMSE = %.6f\n', ...
                pump,tank,step_sizes(i)*100,RMSE);

        end
    end
end
%they fit well

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



function y = transferStepResponse(params,t,tank,pump)

    Kij = params(1);
    Tmain = params(2);

    if tank == pump
        % Diagonal elements: first-order transfer function
        G = Kij / tf([Tmain 1], 1);
    else
        % Off-diagonal elements: second-order transfer function
        Textra = params(3);

        G = Kij / ...
            (tf([Tmain 1], 1) * tf([Textra 1], 1));
    end

    y = squeeze(step(G,t));

end
