function [x, info] = qpsolver(H, g, l, u, A, bl, bu, xinit)
    % Convert bl <= A*x <= bu into inequality form:
    % [ A; -A ] * x <= [ bu; -bl ]

    % If A is empty, then there are no inequality constraints.
    if isempty(A)
        Aineq = [];
        bineq = [];
    else
        Aineq = [A; -A];
        bineq = [bu; -bl];
    end

    % If xinit is not provided, set it to an empty array.
    if nargin < 8 || isempty(xinit)
        xinit = [];
    end

    % Internal quadprog solver options.
    optimalityTolerance = 1e-9;
    maximumIterations = 500;
    solverDisplay = "final"; % Use "iter" to display solver progress, or "off" to suppress output, or "final" to display only the final output.

    options = optimoptions("quadprog", ...
        "OptimalityTolerance", optimalityTolerance, ...
        "MaxIterations", maximumIterations, ...
        "Display", solverDisplay);

    [x, fval, exitflag] = quadprog( ...
        H, g, Aineq, bineq, [], [], l, u, xinit, options);

    info.exitflag = exitflag;
        % info.exitflag > 0  success
        % info.exitflag == 0 stopped without clear success
        % info.exitflag < 0 infeasible or solver failure
    info.fval = fval;
end