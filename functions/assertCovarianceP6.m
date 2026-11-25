function assertCovarianceP6(Q, name)
    name = string(name);
    assert(size(Q, 1) == size(Q, 2), name + " must be square.");
    scale = max(1, norm(Q, 2));
    assert(isreal(Q) && all(isfinite(Q(:))), name + " must be real and finite.");
    assert(norm(Q - Q', "fro") <= 1e-10 * scale, name + " must be symmetric.");
    minimumEigenvalue = min(eig((Q + Q') / 2));
    assert(minimumEigenvalue >= -1e-9 * scale, name + " must be positive semidefinite.");
end
