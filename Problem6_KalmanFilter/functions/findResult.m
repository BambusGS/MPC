function r = findResult(results, scenario, plantType, modelName)
    % Select exactly one simulation result; error if there is none or more than one.
    % (A missing result would otherwise expand to zero function arguments.)
    mask = strcmp({results.scenario}, char(scenario)) ...
         & strcmp({results.plant}, char(plantType)) ...
         & strcmp({results.modelName}, char(modelName));
    if nnz(mask) ~= 1
        error('findResult:notUnique', 'Found %d results for %s / %s / %s (expected 1).', ...
            nnz(mask), char(scenario), char(plantType), char(modelName));
    end
    r = results(mask);
end
