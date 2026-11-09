function r = register_rigid(source, target, varargin)
%REGISTER_RIGID  Rigid registration of SOURCE onto TARGET by ICP.
%   Options: 'method' ('point_to_point' or 'point_to_plane'),
%   'max_iterations' (50), 'tolerance' (1e-10), 'max_distance' (0 keeps
%   every pair), 'initial' (a starting 4x4 matrix, [] for identity, or
%   'auto' for the best of the 24 principal-frame candidates),
%   'symmetric_check' (false), 'source_region' and 'target_region' ('' for
%   all) and 'adaptive_trim' (0: drop pairs beyond this multiple of the
%   median pair distance). Returns a struct: matrix (4x4, apply with
%   mio.transform), initial_rms, rms, iterations, num_pairs, converged,
%   num_rejected, chosen_candidate (0-based, -1 without 'auto'),
%   candidate_rms and warnings.

    o = mio.internal.opts(varargin, {'method', 'point_to_point', 'max_iterations', 50, ...
        'tolerance', 1e-10, 'max_distance', 0, 'initial', [], 'symmetric_check', false, ...
        'source_region', '', 'target_region', '', 'adaptive_trim', 0}, 'mio.register_rigid');
    mode = 'matrix';
    if ischar(o.initial)
        mode = o.initial;
        initial = [];
    elseif isempty(o.initial)
        initial = [];
    else
        initial = double(o.initial);
    end
    r = meshioplusplus_mex('mio_register_rigid', source, target, ...
        mio.internal.str(o.method, 'method'), double(o.max_iterations), double(o.tolerance), ...
        double(o.max_distance), initial, mode, logical(o.symmetric_check), ...
        mio.internal.str(o.source_region, 'source_region'), ...
        mio.internal.str(o.target_region, 'target_region'), double(o.adaptive_trim));
    if isempty(r.warnings)
        r.warnings = cell(0, 1);
    else
        r.warnings = strsplit(r.warnings, char(10));
        r.warnings = r.warnings(:);
    end
end
