function [Xpred, Ypred] = mpc_predict(A, B, C, G, xhat, what, U)
%MPC_PREDICT Recursively predict states and outputs over the MPC horizon

% Infer the horizon from the stacked input sequence and check its dimensions
n = size(A, 1);
m = size(B, 2);
p = size(C, 1);
N = numel(U) / m;
xhat = xhat(:);
what = what(:);
U = U(:);
if size(A, 2) ~= n || size(B, 1) ~= n || size(C, 2) ~= n || size(G, 1) ~= n || N < 1 || N ~= floor(N) || numel(xhat) ~= n || numel(what) ~= size(G, 2)
    error("mpc_predict:inputDimensions", "Prediction inputs have incompatible dimensions.");
end

% Include the current disturbance estimate only in the first step, then propagate nominal dynamics as specified by the MPC model
Xpred = zeros(n, N);
Ypred = zeros(p, N);
x = A * xhat + B * U(1:m) + G * what; % The first predicted state is the free response to the current state, disturbance, and first planned input

for i = 1:N
    if i > 1
        columns = (i - 1) * m + (1:m);
        x = A * x + B * U(columns); % Next predicted state is noise-free (best guess is 0-mean disturbance)
    end
    Xpred(:, i) = x;
    Ypred(:, i) = C * x;
end

end