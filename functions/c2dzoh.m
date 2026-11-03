function [A_d, B_d] = c2dzoh(A, B, Ts)
    % Discretize a continuous-time state-space model
    % using zero-order hold.
    %
    % Continuous-time model:
    %
    %     x_dot = A*x + B*u
    %
    % Discrete-time model:
    %
    %     x(k+1) = A_d*x(k) + B_d*u(k)
    %
    % Inputs:
    %   A: Continuous-time state matrix
    %   B: Continuous-time input matrix
    %   Ts: Sampling time [s]
    %
    % Outputs:
    %   A_d: Discrete-time state matrix
    %   B_d: Discrete-time input matrix

    % Get the number of states and inputs
    [numberOfStates, numberOfInputs] = size(B);

    % Construct the augmented matrix
    augmentedMatrix = [
        A, B;
        zeros(numberOfInputs, numberOfStates + numberOfInputs)
    ];

    % Calculate the matrix exponential
    augmentedExponential = expm(augmentedMatrix * Ts);

    % Extract the discrete-time matrices
    A_d = augmentedExponential( ...
        1:numberOfStates, ...
        1:numberOfStates);

    B_d = augmentedExponential( ...
        1:numberOfStates, ...
        numberOfStates + 1:end);
end