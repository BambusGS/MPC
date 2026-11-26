function [Ustar, uApply, info] = mpc_compute(ctrl, xhat, what, R, uPrev)
%MPC_COMPUTE Solve the unconstrained finite-horizon MPC problem

% Normalize online signals to column vectors and reject mismatched inputs
xhat = xhat(:);
what = what(:);
R = R(:);
uPrev = uPrev(:);
if numel(xhat) ~= size(ctrl.Phi_x, 2) || numel(what) ~= ctrl.nDisturbances || numel(R) ~= size(ctrl.Gamma, 1) || numel(uPrev) ~= ctrl.nInputs
    error("mpc_compute:inputDimensions", "Online MPC inputs have incompatible dimensions.");
end

% Form the online gradient using the free response and previous input
b = ctrl.Phi_x * xhat + ctrl.Phi_w * what;
f = -(ctrl.Gamma' * ctrl.Qbar * (R - b)) - (ctrl.Lambda' * ctrl.Sbar * ctrl.I0 * uPrev);

% Solve H*U = -f directly; this step is valid for the unconstrained problem only though
Ustar = ctrl.H \ (-f);
uApply = Ustar(1:ctrl.nInputs);
info.b = b;
info.f = f;
end
