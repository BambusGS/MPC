    function lin = computeJacobians(op, p)
    % Continuous-time linearization of the four-tank system at a steady state
    %
    % First-order Taylor expansion of the model around the steady state op:
    %   dX/dt = A*X + B*U + E*D,   Y = C*X,   Z = Cz*X
    % with the deviation variables X = x - xs, U = u - us, D = d - ds.
    %
    % Inputs:
    %   op: Steady state created by computeSteadyState()
    %   p: Parameter structure created by createParameters()
    %
    % Output:
    %   lin: Structure with the fields
    %       A: df/dx [1/s]
    %       B: df/du [g/cm^3]
    %       E: df/dd [g/cm^3]
    %       C: dg/dx, measured heights [1/cm^2]
    %       Cz: dh/dx, controlled heights [1/cm^2]
    %       T: Time constants of the four tanks [s]
    %       op: The steady state the model belongs to

    a = p.a;         % Pipe areas [cm^2]
    At = p.A;        % Tank areas [cm^2]
    gamma = p.gamma; % Valve positions [-]
    g = p.g;         % Gravity [cm/s^2]
    rho = p.rho;     % Water density [g/cm^3]

    % Time constants T_i = (A_i/a_i)*sqrt(2*h_i/g) [s]
    T = (At ./ a) .* sqrt(2 * op.h / g);

    % df/dx: each tank loses its own outflow, tanks 3 and 4 feed tanks 1 and 2
    lin.A = [-1/T(1),  0,       1/T(3),  0;
              0,      -1/T(2),  0,       1/T(4);
              0,       0,      -1/T(3),  0;
              0,       0,       0,      -1/T(4)];

    % df/du and df/dd: the flows enter the mass balances linearly
    lin.B = rho * [gamma(1),     0;
                   0,            gamma(2);
                   0,            1 - gamma(2);
                   1 - gamma(1), 0];
    lin.E = rho * [0, 0;
                   0, 0;
                   1, 0;
                   0, 1];

    % dg/dx and dh/dx: only the heights of tanks 1 and 2 are measured and controlled
    lin.C = [1/(rho*At(1)), 0,             0, 0;
             0,             1/(rho*At(2)), 0, 0];
    lin.Cz = lin.C;

    lin.T = T;
    lin.op = op;
end
