function h = graphicsArray(n)
% GRAPHICSARRAY  Preallocated array of graphics handles.
    if exist('gobjects', 'file') || exist('gobjects', 'builtin')
        h = gobjects(1, n);
    else
        h = zeros(1, n);
    end
end
