function H = computeMarkovParameters(A, B, C, N)
    % Markov parameters (impulse response coefficients) of a discrete-time model
    %
    %   H_0 = D = 0,   H_i = C*A^(i-1)*B,   i = 1, ..., N
    %
    % Inputs:
    %   A: Discrete-time state matrix
    %   B: Discrete-time input matrix
    %   C: Output matrix
    %   N: Number of Markov parameters after H_0
    %
    % Output:
    %   H: ny-by-nu-by-(N+1) array, H(:, :, i+1) = H_i

    H = zeros(size(C, 1), size(B, 2), N + 1);

    % Propagate A^(i-1)*B instead of computing matrix powers
    AiB = B;
    for i = 1:N
        H(:, :, i + 1) = C * AiB;
        AiB = A * AiB;
    end
end
