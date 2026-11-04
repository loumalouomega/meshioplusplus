function r = cut_quadrature(mesh, varargin)
%CUT_QUADRATURE  Quadrature rules on the part of a mesh that level sets select.
%   R = mio.cut_quadrature(MESH, ...) builds cut-cell quadrature rules through
%   CutCells: rules of degree 'order' on the part mio.cut_level_set would
%   give, as a point cloud with point data 'quadrature:weight',
%   'quadrature:reference', 'quadrature:rule' and 'quadrature:parent_cell'.
%   Options (as mio.cut_level_set, plus 'order'):
%     'level_sets' ({}), 'planes' ([]), 'selection' (''),
%     'order'      the polynomial degree the rules integrate exactly (2),
%     'include_uncut' (true), 'backend' ('auto'), 'triangulate' (false),
%     'record_parent_ids' (true)
%   R is a struct: mesh, num_rules, num_points, sum_weights, dimension and
%   backend.

    o = mio.internal.opts(varargin, {'level_sets', {}, 'planes', [], 'selection', '', ...
        'order', 2, 'include_uncut', true, 'backend', 'auto', 'triangulate', false, ...
        'record_parent_ids', true}, 'mio.cut_quadrature');
    r = meshioplusplus_mex('mio_cut_quadrature', mesh, ...
        mio.internal.cellstr_arg(o.level_sets, 'level_sets'), cut_planes(o.planes), ...
        mio.internal.str(o.selection, 'selection'), double(o.order), ...
        logical(o.include_uncut), mio.internal.str(o.backend, 'backend'), ...
        logical(o.triangulate), logical(o.record_parent_ids));
end

function flat = cut_planes(planes)
%   Six numbers per plane, flat: [origin1 normal1 origin2 normal2 ...].
    if isempty(planes)
        flat = zeros(0, 1);
        return;
    end
    if iscell(planes)
        parts = cellfun(@(p) double(p(:)), planes(:), 'UniformOutput', false);
        flat = vertcat(parts{:});
    elseif isvector(planes)
        flat = double(planes(:));
    elseif size(planes, 2) == 6
        flat = reshape(double(planes).', [], 1);
    else
        flat = double(planes(:));
    end
    if mod(numel(flat), 6) ~= 0
        error('meshioplusplus:argument', ...
              'meshio++: cut: planes must hold six numbers each (origin, then normal)');
    end
end
