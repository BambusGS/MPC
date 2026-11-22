function ident = sharedTimeConstants(estimatedParams, stepIndex)
    % Identified model of Problem 4 with the time constants shared as in the lecture
    %
    %   G11 = K11/(T1 s + 1)                G12 = K12/((T3 s + 1)(T1 s + 1))
    %   G21 = K21/((T4 s + 1)(T2 s + 1))    G22 = K22/(T2 s + 1)
    %
    % problem4_main.m fits T1 twice (in G11 and G12) and T2 twice (in G21 and G22).
    % In the off-diagonal elements the order of Tmain and Textra is arbitrary, so
    % the lag closest to the diagonal time constant is taken as T1 (T2), and the
    % other lag as T3 (T4). The shared T1 (T2) is the mean of the two fits.
    %
    % Inputs:
    %   estimatedParams: estimated_params from problem4_main.m,
    %       estimatedParams(tank, pump, :, step) = [K, Tmain, Textra]
    %   stepIndex: 1, 2 or 3 for the 10 %, 25 % or 50 % fit
    %
    % Output:
    %   ident: Structure with the fields
    %       K: 2x2 steady-state gains [cm/(cm^3/s)]
    %       T: [T1; T2; T3; T4] [s]

    P = estimatedParams(:, :, :, stepIndex);

    ident.K = P(:, :, 1);

    T11 = P(1, 1, 2);
    T22 = P(2, 2, 2);
    lags12 = squeeze(P(1, 2, 2:3));
    lags21 = squeeze(P(2, 1, 2:3));

    [~, i1] = min(abs(lags12 - T11));
    [~, i2] = min(abs(lags21 - T22));

    T1 = mean([T11, lags12(i1)]);
    T2 = mean([T22, lags21(i2)]);
    T3 = lags12(3 - i1);
    T4 = lags21(3 - i2);

    ident.T = [T1; T2; T3; T4];
end
