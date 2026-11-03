function y = systemSense(x, p)
    % Return the two measured bottom-tank heights
    %
    % Inputs:
    %   x: State vector containing tank masses [g]
    %   p: Parameter structure
    %
    % Output:
    %   y: Measured tank 1 and tank 2 heights [cm]

    % Unpack parameters
    A = p.A;       % Tank areas [cm^2]
    rho = p.rho;   % Water density [g/cm^3]

    % Only tanks 1 and 2 have sensors
    y = x(1:2) ./ (rho .* A(1:2));
end