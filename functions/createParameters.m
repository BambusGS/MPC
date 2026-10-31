function p = createParameters()
    % Four-tank process parameters

    p.a = [1.2272; 1.2272; 1.2272; 1.2272];       % Pipe areas [cm^2]
    p.A = [380.1327; 380.1327; 380.1327; 380.1327]; % Tank areas [cm^2]

    p.gamma = [0.60; 0.70]; % Valve positions [-]

    p.g = 981;   % Gravity [cm/s^2]
    p.rho = 1.0; % Water density [g/cm^3]

    % Disturbance model parameters
    p.FMean = [50; 100];       % Mean F3/F4 disturbance flows [cm^3/s]
    p.piecewiseSigma = [10; 20]; % Piecewise disturbance standard deviations [cm^3/s]
    p.disturbanceRate = [0.05; 0.05]; % Mean-reversion rates [1/s most probably]
    p.disturbanceSigma = [5; 5];    % Diffusion coefficients
    p.measurementVariance = 1; % Variance of each level sensor [cm^2]
    p.measurementCovariance = kron(p.measurementVariance, eye(2)); % y1/y2 sensor noise covariance [cm^2]
end