function detectable = isDetectable(A, C)
    % PBH test: (C, A) is detectable if rank([lambda*I - A; C]) = n for every
    % eigenvalue lambda of A with |lambda| >= 1 (including the integrators).

    n = size(A, 1);
    lambda = eig(A);
    detectable = true;
    for i = 1:n
        if abs(lambda(i)) >= 1 - 1e-9 && rank([lambda(i) * eye(n) - A; C], 1e-9) < n
            detectable = false;
        end
    end
end
