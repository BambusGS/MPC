function sim = createSimSettings(p)
    % Simulation and controller settings

    % Simulation time settings
    sim.t0 = 0;       % Start time [s]
    sim.tf = 20 * 60; % End time [s]
    sim.Ts = 10;      % Sampling time [s] used for controller and measurements
    sim.dt = 0.02;     % Time step [s] used for numerical integration of the SDE

    % Time vector and number of samples
    sim.t = (sim.t0:sim.Ts:sim.tf).';
    sim.N = numel(sim.t);

    % Initial tank masses [g]
    sim.initialState = zeros(4, 1);
    sim.initialDisturbance = p.FMean;

    % Pump flow rates [cm^3/s]
    sim.F1 = 300;
    sim.F2 = 300;

    % Additional flow parameters, if needed later
    sim.F3_mean = p.FMean(1);
    sim.F4_mean = p.FMean(2);

    % Constant input matrix for open-loop simulation
    sim.input = repmat([sim.F1; sim.F2], 1, sim.N);

    % Reference heights for tanks 1 and 2 [cm]
    sim.reference = [30; 30];

    % PI controller parameters
    sim.Kp = 0.2;
    sim.Tinteg = 5000;

    % Initial integral term
    sim.integralTerm = zeros(2, 1);

    sim.randomSeed = 68;
end