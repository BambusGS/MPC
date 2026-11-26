function ctrl = mpc_design(A, B, C, G, N, Q, S)
%MPC_DESIGN Build prediction matrices and gains for unconstrained MPC
% Inputs use y = C*x and x(k+1) = A*x + B*u + G*w
% N is the prediction horizon in steps, Q is the output weight, and S is the input weight

% Check dimensions once during design
n = size(A, 1);
m = size(B, 2);
p = size(C, 1);
q = size(G, 2);
if size(A, 2) ~= n || size(B, 1) ~= n || size(C, 2) ~= n || size(G, 1) ~= n
    error("mpc_design:modelDimensions", "A, B, C, and G have incompatible dimensions.");
end
if ~isscalar(N) || N < 1 || N ~= floor(N)
    error("mpc_design:horizon", "N must be a positive integer.");
end
if ~isequal(size(Q), [p, p]) || ~isequal(size(S), [m, m])
    error("mpc_design:weightDimensions", "Q and S must match output and input dimensions.");
end
if norm(Q - Q', "fro") > 1e-10 || norm(S - S', "fro") > 1e-10
    error("mpc_design:weightSymmetry", "Q and S must be symmetric.");
end


% Stack the effects of the initial state, current disturbance, and planned inputs on the N future outputs
Phi_x = zeros(N * p, n);
Phi_w = zeros(N * p, q);
Gamma = zeros(N * p, N * m);
for predStep = 1:N
    outputRows = (predStep - 1) * p + (1:p);   % Select the p rows corresponding to the current predicted output
    Phi_x(outputRows, :) = C * A^predStep;
    Phi_w(outputRows, :) = C * A^(predStep - 1) * G;
    for inputStep = 1:predStep
        inputColumns = (inputStep - 1) * m + (1:m);
        Gamma(outputRows, inputColumns) = C * A^(predStep - inputStep) * B;
    end
end


% Lambda converts the planned inputs to input increments, including the difference between the first planned input and the previous applied input
Lambda = kron(eye(N), eye(m)); % The [-I I 0 0 ... matrix row]
for predStep = 2:N
    rows = (predStep - 1) * m + (1:m);
    columns = (predStep - 2) * m + (1:m);
    Lambda(rows, columns) = -eye(m);
end
I0 = [eye(m); zeros((N - 1) * m, m)]; % I0 selects the first planned input from the N stacked inputs, so that u(k) = I0 * u_bar(k)
Qbar = kron(eye(N), Q); % Makes block-diag matrix
Sbar = kron(eye(N), S); % Makes block-diag matrix


% The Hessian is fixed for a given model and tuning, so calculate it once
H = Gamma' * Qbar * Gamma + Lambda' * Sbar * Lambda;
H = (H + H') / 2;
[~, cholFlag] = chol(H); % Use Cholesky factorization to check that the Hessian is positive definite
if cholFlag ~= 0
    error("mpc_design:hessianNotPositiveDefinite", "The MPC Hessian must be positive definite.");
end

% Store first-input gains for the equivalent affine control law
ctrl.A = A;
ctrl.B = B;
ctrl.C = C;
ctrl.G = G;
ctrl.N = N;
ctrl.nInputs = m;
ctrl.nOutputs = p;
ctrl.nDisturbances = q;
ctrl.Gamma = Gamma;
ctrl.Phi_x = Phi_x;
ctrl.Phi_w = Phi_w;
ctrl.Lambda = Lambda;
ctrl.I0 = I0;
ctrl.Qbar = Qbar;
ctrl.Sbar = Sbar;
ctrl.H = H;

% Gain matrices extraction
selector = [eye(m), zeros(m, (N - 1) * m)];
ctrl.Lx = -selector * (H \ (Gamma' * Qbar * Phi_x));
ctrl.Lw = -selector * (H \ (Gamma' * Qbar * Phi_w));
ctrl.LR = selector * (H \ (Gamma' * Qbar));
ctrl.Lu = selector * (H \ (Lambda' * Sbar * I0));

end
