function r = clean(mesh, varargin)
%CLEAN  Weld points, drop orphans, degenerate and duplicate cells.
%   R = mio.clean(MESH, ...) with options 'weld' (false: fuse coincident
%   points within 'atol'), 'atol' (1e-12), 'remove_orphans' (true),
%   'drop_degenerate' (true) and 'drop_duplicate_cells' (true). Returns a
%   struct: mesh, points_welded, points_removed_orphan,
%   cells_dropped_degenerate and cells_dropped_duplicate.

    o = mio.internal.opts(varargin, {'weld', false, 'atol', 1e-12, 'remove_orphans', true, ...
        'drop_degenerate', true, 'drop_duplicate_cells', true}, 'mio.clean');
    r = meshioplusplus_mex('mio_clean', mesh, logical(o.weld), double(o.atol), ...
        logical(o.remove_orphans), logical(o.drop_degenerate), logical(o.drop_duplicate_cells));
end
