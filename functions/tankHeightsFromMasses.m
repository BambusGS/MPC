function heights = tankHeightsFromMasses(masses, p)
    % Convert four tank masses or mass deviations to heights [cm].
    %
    % masses must have one row per tank and may contain multiple time samples.
    % The returned heights have the same shape as masses.

    tankAreas = p.A(:);
    if size(masses, 1) ~= numel(tankAreas)
        error("tankHeightsFromMasses:InvalidMassShape", ...
            "Mass data must have one row for each of the four tanks.");
    end

    heights = masses ./ (p.rho .* tankAreas);
end
