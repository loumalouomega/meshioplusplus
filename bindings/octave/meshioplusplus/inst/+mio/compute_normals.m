function r = compute_normals(mesh, varargin)
%COMPUTE_NORMALS  Point and cell normals of a surface, optionally split at creases.
%   Writes 'normals' point data (n x 3) and, with 'cell_normals', cell data.
%   Options: 'point_normals' (true), 'cell_normals' (false), 'weight'
%   ('angle' or 'area'), 'split_angle' ([] for one smooth normal per point,
%   or the largest dihedral angle in degrees still smooth; split copies are
%   appended after the original points), 'record_parent_ids' (false: attach
%   'normals:parent_point') and 'region' (''). Returns a struct: mesh,
%   quality, num_isolated, num_undefined, num_degenerate, num_split_points
%   and num_added_points.

    o = mio.internal.opts(varargin, {'point_normals', true, 'cell_normals', false, ...
        'weight', 'angle', 'split_angle', [], 'record_parent_ids', false, 'region', ''}, ...
        'mio.compute_normals');
    if isempty(o.split_angle)
        split_angle = [];
    else
        split_angle = double(o.split_angle);
    end
    r = meshioplusplus_mex('mio_compute_normals', mesh, logical(o.point_normals), ...
        logical(o.cell_normals), mio.internal.str(o.weight, 'weight'), split_angle, ...
        logical(o.record_parent_ids), mio.internal.str(o.region, 'region'));
end
