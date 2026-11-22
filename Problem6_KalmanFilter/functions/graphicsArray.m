function h = graphicsArray(n)
    % Preallocated array for graphics handles (gobjects in MATLAB, zeros in Octave).
    % MATLAB cannot store graphics objects in a double array.
    if exist('gobjects', 'file') || exist('gobjects', 'builtin')
        h = gobjects(1, n);
    else
        h = zeros(1, n);
    end
end
