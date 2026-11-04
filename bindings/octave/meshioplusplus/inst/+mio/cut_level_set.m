function r = cut_level_set(mesh, varargin)
%CUT_LEVEL_SET  Cut a mesh by level sets.
%   R = mio.cut_level_set(MESH, ...) returns the part of MESH that level sets
%   select, through CutCells (doc/cut.md): the pieces of the cells they cross
%   on one side of them or on their zero sets, and with 'include_uncut' the
%   whole cells inside. Point data is interpolated; cell data and cell
%   regions follow each piece's parent ('cut:parent_cell'). A library without
%   CutCells fails naming MESHIOPLUSPLUS_WITH_CUTCELLS. Options:
%     'level_sets'        cellstr of scalar point-data arrays that are level sets ({})
%     'planes'            plane level sets, each [origin, normal] (six
%                         numbers): an n x 6 matrix (one plane per row) or a
%                         cell array of 6-vectors; named plane0, plane1, ... ([])
%     'selection'         the part, such as 'phi < 0', 'phi = 0' or
%                         'a < 0 and b > 0'; empty is '<first> < 0' ('')
%     'include_uncut'     keep the whole cells inside the part (true)
%     'backend'           'auto', 'lut' or 'quadrays' ('auto')
%     'degree'            quadrays only: the pieces' Lagrange degree (2)
%     'triangulate'       lut only: split the pieces into simplices (false)
%     'record_parent_ids' write 'cut:parent_cell' (true)
%   R is a struct: mesh, num_cells, num_cut_cells, num_uncut_cells,
%   num_pieces, dimension, measure and backend.

    o = mio.internal.opts(varargin, {'level_sets', {}, 'planes', [], 'selection', '', ...
        'include_uncut', true, 'backend', 'auto', 'degree', 2, 'triangulate', false, ...
        'record_parent_ids', true}, 'mio.cut_level_set');
    r = meshioplusplus_mex('mio_cut_level_set', mesh, ...
        mio.internal.cellstr_arg(o.level_sets, 'level_sets'), cut_planes(o.planes), ...
        mio.internal.str(o.selection, 'selection'), logical(o.include_uncut), ...
        mio.internal.str(o.backend, 'backend'), double(o.degree), logical(o.triangulate), ...
        logical(o.record_parent_ids));
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
