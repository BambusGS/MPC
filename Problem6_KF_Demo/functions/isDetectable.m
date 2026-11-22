function detectable = isDetectable(A, C)
% ISDETECTABLE  PBH detectability test of the pair (C, A).
    n = size(A, 1);
    lambda = eig(A);
    detectable = true;
    for i = 1:n
        if abs(lambda(i)) >= 1 - 1e-9 && rank([lambda(i) * eye(n) - A; C], 1e-9) < n
            detectable = false;
        end
    end
end
