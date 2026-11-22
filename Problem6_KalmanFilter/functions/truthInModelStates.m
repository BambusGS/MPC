function X = truthInModelStates(p, op, plant, model)
    % True plant trajectory expressed in the states of a Kalman filter model,
    % in the physical units used for plotting (model.scale)
    %
    % Problem 5 model: [h1..h4; F3; F4 (; N3; N4)] deviations, N = FMean - d_s.
    % Problem 4 model: x1, x2 = h1, h2 deviations; x3 = dq3/(1 - gamma2) and
    %   x4 = dq4/(1 - gamma1), the pump-flow equivalents of the upper-tank
    %   outflows; N = (FMean - d_s) in the same units. These equivalents use
    %   gamma and are for plotting only; the filter does not know them.

    g1 = p.gamma(1);
    g2 = p.gamma(2);
    dh = plant.h - op.h;
    dq = plant.q - op.q;
    dF = plant.F - op.d;
    dN = plant.FMean - op.d;

    if model.problem == 5
        X = [dh; dF];
        if model.offsetFree
            X = [X; dN];
        end
    else
        X = [dh(1:2, :); dq(3, :) / (1 - g2); dq(4, :) / (1 - g1)];
        if model.offsetFree
            X = [X; dN(1, :) / (1 - g2); dN(2, :) / (1 - g1)];
        end
    end
end
