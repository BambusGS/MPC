function [Ustar, uApply, info] = mpc_compute_qp(ctrl, xhat, what, R, uPrev, constraints)
%MPC_COMPUTE_QP Solve the same MPC cost with linear constraints or bounds constraints may contain l, u, A, bl, bu, and xinit for qpsolver

% Build the same online gradient as the unconstrained controller
xhat = xhat(:);
what = what(:);
R = R(:);
uPrev = uPrev(:);
if numel(xhat) ~= size(ctrl.Phi_x, 2) || numel(what) ~= ctrl.nDisturbances || numel(R) ~= size(ctrl.Gamma, 1) || numel(uPrev) ~= ctrl.nInputs
    error("mpc_compute_qp:inputDimensions", "Online MPC inputs have incompatible dimensions.");
end

b = ctrl.Phi_x * xhat + ctrl.Phi_w * what;
f = -(ctrl.Gamma' * ctrl.Qbar * (R - b)) - (ctrl.Lambda' * ctrl.Sbar * ctrl.I0 * uPrev);

% Convert optional constraint fields to the qpsolver's expected arguments
if nargin < 6 || isempty(constraints)
    constraints = struct();
end
l = [];
u = [];
A = [];
bl = [];
bu = [];
xinit = [];
if isfield(constraints, "l"), l = constraints.l; end
if isfield(constraints, "u"), u = constraints.u; end
if isfield(constraints, "A"), A = constraints.A; end
if isfield(constraints, "xinit"), xinit = constraints.xinit; end
if ~isempty(A)
    if ~isfield(constraints, "bl") || ~isfield(constraints, "bu")
        error("mpc_compute_qp:constraintBounds", "constraints.bl and constraints.bu are required when constraints.A is nonempty.");
    end
    bl = constraints.bl;
    bu = constraints.bu;
end

% qpsolver handles bounds and two-sided linear constraints on the stacked U
[Ustar, info] = qpsolver(ctrl.H, f, l, u, A, bl, bu, xinit);
if info.exitflag <= 0
    error("mpc_compute_qp:solverFailure", "QP solver failed with exit flag %d.", info.exitflag);
end
uApply = Ustar(1:ctrl.nInputs);
info.b = b;
info.f = f;

end