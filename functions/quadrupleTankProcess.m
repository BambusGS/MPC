function xdot = quadrupleTankProcess(~, x, u, d, p)
    % Four tank system model: dx/dt = f(t,x,u,p)
    %
    % This function implements a differential equation model
    % for the four-tank system.
    %
    % Inputs:
    %   t: Current simulation time [s]
    %   x: State vector containing the mass in each tank [g]
    %   u: Input vector containing the pump flow rates [cm^3/s]
    %   d: Disturbance vector containing F3 and F4 [cm^3/s]
    %   p: Parameter structure created by createParameters()
    %
    % Output:
    %   xdot: Time derivative of the state vector [g/s]

    % Unpack states and inputs
    m = x; % Mass of liquid in each tank [g]
    F = u;         % Flow rates from the pumps [cm^3/s]
    disturbance = d;

    % Unpack parameters from the parameter structure
    a = p.a; % Pipe cross-sectional areas [cm^2]
    A = p.A; % Tank cross-sectional areas [cm^2]
    gamma = p.gamma; % Valve positions [-]
    g = p.g; % Acceleration due to gravity [cm/s^2]
    rho = p.rho; % Density of water [g/cm^3]

    
    % Inflows to the four tanks [cm^3/s]
    qin = zeros(4, 1); 
    qin(1, 1) = gamma(1) * F(1); 
    qin(2, 1) = gamma(2) * F(2); 
    qin(3, 1) = (1 - gamma(2)) * F(2); 
    qin(4, 1) = (1 - gamma(1)) * F(1); 
    
    % Outflows from the tanks [cm^3/s]
    h = m ./ (rho .* A); % Liquid level in each tank [cm]
    qout = a .* sqrt(2 * g * h); % Outflow from each tank [cm3/s]
    % qout = zeros(4, 1);
    % qout(1, 1) = a(1) * sqrt(2 * g * h(1)); 
    % qout(2, 1) = a(2) * sqrt(2 * g * h(2)); 
    % qout(3, 1) = a(3) * sqrt(2 * g * h(3)); 
    % qout(4, 1) = a(4) * sqrt(2 * g * h(4)); 
    
    % Differential equations for the mass balances [g/s]
    xdot = zeros(4, 1); 
    xdot(1, 1) = rho * (qin(1, 1) + qout(3, 1) - qout(1, 1));
    xdot(2, 1) = rho * (qin(2, 1) + qout(4, 1) - qout(2, 1));
    xdot(3, 1) = rho * (qin(3, 1) - qout(3, 1) + disturbance(1));
    xdot(4, 1) = rho * (qin(4, 1) - qout(4, 1) + disturbance(2));
end