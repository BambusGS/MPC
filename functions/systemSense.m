function y = systemSense(x, p)
    % Return the two measured bottom-tank heights
    %
    % Inputs:
    %   x: State vector containing tank masses [g]
    %   p: Parameter structure
    %
    % Output:
    %   y: Measured tank 1 and tank 2 heights [cm]

    % Only tanks 1 and 2 have sensors
    heights = tankHeightsFromMasses(x, p);
    y = heights(1:2);
end