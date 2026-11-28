function r = subsample_points(mesh, count, varargin)
%SUBSAMPLE_POINTS  Subsample a point cloud to an exact count budget.
%   R = mio.subsample_points(MESH, COUNT, ...) selects COUNT points using
%   'method' ('farthest', 'grid', or 'random'), 'seed' (0), 'start' ([] for
%   seeded draw, or 1-based index), 'bounds' ([] or 1x6 [xlo ylo zlo xhi yhi zhi]),
%   'max_cells' (16000000), and 'record_ids' (false).
%   Returns a struct: mesh, point_map (n x 1, 1-based), and warnings.

    o = mio.internal.opts(varargin, {'method', 'farthest', 'seed', 0, ...
        'start', [], 'bounds', [], 'max_cells', 16000000, 'record_ids', false}, ...
        'mio.subsample_points');
    s = o.start;
    if ~isempty(s)
        s = double(s);
    end
    b = o.bounds;
    if ~isempty(b)
        b = double(b(:));
    end
    r = meshioplusplus_mex('mio_subsample_points', mesh, double(count), ...
        mio.internal.str(o.method, 'method'), double(o.seed), s, b, ...
        double(o.max_cells), logical(o.record_ids));
    if isempty(r.warnings)
        r.warnings = {};
    else
        r.warnings = strsplit(r.warnings, char(10));
    end
end
