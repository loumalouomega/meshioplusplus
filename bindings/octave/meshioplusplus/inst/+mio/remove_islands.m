function r = remove_islands(mesh, varargin)
%REMOVE_ISLANDS  Remove the small connected components of a mesh.
%   R = mio.remove_islands(MESH, ...) drops a scan's crumbs (roadmap 7.8.1):
%   components joined by shared edges ('connectivity', 'edge') or points
%   ('point') whose bounding-box diagonal ('min_diagonal', 0.01), area
%   ('min_area', 0) or cell count ('min_cells', 0) is below the threshold,
%   'relative' (true) to the whole mesh. 'dry_run' (false) labels the cells
%   instead. R is a struct: mesh, num_components, num_removed,
%   num_cells_removed and num_points_removed.

    o = mio.internal.opts(varargin, {'connectivity', 'edge', 'min_diagonal', 0.01, ...
        'min_area', 0, 'min_cells', 0, 'relative', true, 'dry_run', false}, ...
        'mio.remove_islands');
    conn = find(strcmp({'edge', 'point'}, o.connectivity));
    if isempty(conn)
        error('meshioplusplus:argument', ...
            'meshio++: remove_islands: connectivity must be edge or point');
    end
    r = meshioplusplus_mex('mio_remove_islands', mesh, double(conn - 1), ...
        double(o.min_diagonal), double(o.min_area), double(o.min_cells), ...
        logical(o.relative), logical(o.dry_run));
end
