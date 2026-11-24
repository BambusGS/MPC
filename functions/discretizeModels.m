function disc = discretizeModels(lin, aug, p, Ts)
    % Discrete-time models of the three linearized four-tank models
    %
    % Deterministic: X(k+1) = A*X(k) + B*U(k) + E*D(k)
    % Piecewise:     X(k+1) = A*X(k) + B*U(k) + E*w(k),   w(k) ~ N(0, Rdd)
    % SDE:           Xa(k+1) = Aa*Xa(k) + Ba*U(k) + w(k), w(k) ~ N(0, Q)
    % All: Y(k) = C*X(k) + v(k), v(k) ~ N(0, Rvv)
    %
    % Inputs:
    %   lin: Linearized model created by computeJacobians()
    %   aug: Augmented SDE model created by computeAugmentedModel()
    %   p: Parameter structure created by createParameters()
    %   Ts: Sampling time [s]
    %
    % Output:
    %   disc: Structure with the fields
    %       Ts: Sampling time [s]
    %       det: A, B, E, C, Cz (deterministic model)
    %       pw: A, B, E, C, Cz, Rdd, Q = E*Rdd*E' (piecewise constant model)
    %       sde: A, B, C, Cz, Q (SDE model, Q by Van Loan)
    %       Rvv: Measurement noise covariance

    nInputs = size(lin.B, 2);

    % Deterministic model: inputs and disturbances held constant (ZOH)
    [A, BE] = c2dzoh(lin.A, [lin.B, lin.E], Ts);
    disc.det.A = A;
    disc.det.B = BE(:, 1:nInputs);
    disc.det.E = BE(:, nInputs + 1:end);
    disc.det.C = lin.C;
    disc.det.Cz = lin.Cz;

    % Piecewise constant model: the noise is held constant like an input
    disc.pw = disc.det;
    disc.pw.Rdd = diag(p.piecewiseSigma .^ 2);
    disc.pw.Q = disc.pw.E * disc.pw.Rdd * disc.pw.E';

    % SDE model: the noise accumulates over the interval
    [disc.sde.A, disc.sde.B, disc.sde.Q] = c2dsde(aug.A, aug.B, aug.G, Ts);
    disc.sde.C = aug.C;
    disc.sde.Cz = aug.Cz;

    disc.Rvv = p.measurementCovariance;
    disc.Ts = Ts;
end
