function op = computeSteadyState(u, d, p)
    % Steady state of the four-tank system for constant pump and disturbance flows
    %
    % At steady state the outflow of every tank equals its inflow. The flows are
    % solved starting with the top tanks (3 and 4), which do not depend on the
    % others. Torricelli's law q = a*sqrt(2*g*h) is then inverted for the levels.
    %
    % Inputs:
    %   u: Pump flow rates [F1; F2] [cm^3/s]
    %   d: Disturbance flow rates [F3; F4] [cm^3/s]
    %   p: Parameter structure created by createParameters()
    %
    % Output:
    %   op: Structure with the fields
    %       x: Tank masses (state) [g]
    %       y: Measured heights, tanks 1 and 2 [cm]
    %       z: Controlled heights, tanks 1 and 2 [cm]
    %       h: Heights of all four tanks [cm]
    %       q: Outflows of all four tanks [cm^3/s]
    %       u, d: The inputs the steady state belongs to

    a = p.a;         % Pipe areas [cm^2]
    A = p.A;         % Tank areas [cm^2]
    gamma = p.gamma; % Valve positions [-]
    g = p.g;         % Gravity [cm/s^2]
    rho = p.rho;     % Water density [g/cm^3]

    % Outflow = inflow, top tanks first [cm^3/s]
    q = zeros(4, 1);
    q(3) = (1 - gamma(2)) * u(2) + d(1);
    q(4) = (1 - gamma(1)) * u(1) + d(2);
    q(1) = gamma(1) * u(1) + q(3);
    q(2) = gamma(2) * u(2) + q(4);

    % Invert Torricelli's law 
    h = (q ./ a).^2 / (2 * g);
    x = rho .* A .* h; % Masses [g]

    op.x = x;
    op.y = systemSense(x, p);
    op.z = systemOutput(x, p);
    op.h = h;
    op.q = q;
    op.u = u;
    op.d = d;
end
