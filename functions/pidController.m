function [u, integralTerm] = pidController( ...
    integralTerm, reference, measurement, previousInput, Kp, Ti, Ts)

    % PI controller.
    %
    % Inputs:
    %   integralTerm: Integral term from the previous sample
    %   reference: Desired output value
    %   measurement: Measured output value
    %   previousInput: Previous controller input
    %   Kp: Proportional gain
    %   Ti: Integral time constant [s]
    %   Ts: Sampling time [s]
    %
    % Outputs:
    %   u: New controller input
    %   integralTerm: Updated integral term

    % Calculate the control error
    error = reference - measurement;

    % Calculate the new control input
    u = previousInput + Kp * error + integralTerm;

    % Update the integral term
    integralTerm = integralTerm + (Kp * Ts / Ti) * error;
end