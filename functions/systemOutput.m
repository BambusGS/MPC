function z = systemOutput(x, p)
    % Extract controlled variables from the measured system output
    %
    % The controlled variables are the heights of tanks 1 and 2
    %
    % Inputs:
    %   x: State vector containing tank masses [g]
    %   p: Parameter structure
    %
    % Output:
    %   z: Controlled tank heights [cm]

    % The two controlled outputs are the two measured bottom-tank heights
    %  same as the measured outputs y in this case, so we re-use the function
    z = systemSense(x, p);
end