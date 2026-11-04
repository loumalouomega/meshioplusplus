function r = register_rigid(source, target, varargin)
%REGISTER_RIGID  Rigid registration of SOURCE onto TARGET by ICP.
%   Options: 'method' ('point_to_point' or 'point_to_plane'),
%   'max_iterations' (50), 'tolerance' (1e-10), 'max_distance' (0 keeps
%   every pair) and 'initial' (a starting 4x4 matrix, [] for identity).
%   Returns a struct: matrix (4x4, apply with mio.transform), initial_rms,
%   rms, iterations, num_pairs and converged.

    o = mio.internal.opts(varargin, {'method', 'point_to_point', 'max_iterations', 50, ...
        'tolerance', 1e-10, 'max_distance', 0, 'initial', []}, 'mio.register_rigid');
    if isempty(o.initial)
        initial = [];
    else
        initial = double(o.initial);
    end
    r = meshioplusplus_mex('mio_register_rigid', source, target, ...
        mio.internal.str(o.method, 'method'), double(o.max_iterations), double(o.tolerance), ...
        double(o.max_distance), initial);
end
