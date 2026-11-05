function [A_d, B_d, Q_d] = c2dsde(A, B, G, Ts)
    % Discretize a continuous-time linear SDE
    % using zero-order hold for the input and Van Loan's method for the noise.
    %
    % Continuous-time model:
    %
    %     dx = (A*x + B*u)*dt + G*dw,   dw ~ N(0, I*dt)
    %
    % Discrete-time model:
    %
    %     x(k+1) = A_d*x(k) + B_d*u(k) + w(k),   w(k) ~ N(0, Q_d)
    %
    % with Q_d = integral_0^Ts expm(A*t)*G*G'*expm(A'*t) dt.
    %
    % Inputs:
    %   A: Continuous-time state matrix
    %   B: Continuous-time input matrix
    %   G: Continuous-time noise matrix
    %   Ts: Sampling time [s]
    %
    % Outputs:
    %   A_d: Discrete-time state matrix
    %   B_d: Discrete-time input matrix
    %   Q_d: Process noise covariance

    % Deterministic part
    [A_d, B_d] = c2dzoh(A, B, Ts);

    % Construct the Van Loan matrix
    numberOfStates = size(A, 1);
    
    vanLoanMatrix = [
        -A, G * G';
        zeros(numberOfStates), A'
    ];

    % Calculate the matrix exponential
    vanLoanExponential = expm(vanLoanMatrix * Ts);

    % Extract the process noise covariance
    Phi12 = vanLoanExponential(1:numberOfStates, numberOfStates + 1:end);
    Phi22 = vanLoanExponential(numberOfStates + 1:end, numberOfStates + 1:end);
    Q_d = Phi22' * Phi12;
end
