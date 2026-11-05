function aug = computeAugmentedModel(lin, p)
    % Linearized SDE model with the disturbance flows as additional states
    %
    %   dXa = (Aa*Xa + Ba*U)*dt + G*dw,   Y = Ca*Xa,   Z = Cza*Xa
    %   with Xa = [X; D], D = [F3; F4] - FMean.
    %
    % Inputs:
    %   lin: Linearized model created by computeJacobians()
    %   p: Parameter structure created by createParameters()
    %
    % Output:
    %   aug: Structure with the fields
    %       A: [A, E; 0, -diag(a_F)] [1/s]
    %       B: [B; 0]
    %       G: [0; diag(sigma_F)], noise matrix
    %       C: [C, 0], measured heights
    %       Cz: [Cz, 0], controlled heights

    nStates = size(lin.A, 1);
    nDisturbances = size(lin.E, 2);
    nInputs = size(lin.B, 2);

    % Tanks: unchanged, E couples the disturbance states into tanks 3 and 4
    % Disturbances: mean reversion with rate a_F, driven by the Brownian noise
    aug.A = [lin.A, lin.E;
             zeros(nDisturbances, nStates), -diag(p.disturbanceRate)];
    aug.B = [lin.B; zeros(nDisturbances, nInputs)];
    aug.G = [zeros(nStates, nDisturbances); diag(p.disturbanceSigma)];

    % Measurements and outputs depend only on the tank states
    aug.C = [lin.C, zeros(size(lin.C, 1), nDisturbances)];
    aug.Cz = [lin.Cz, zeros(size(lin.Cz, 1), nDisturbances)];
end
