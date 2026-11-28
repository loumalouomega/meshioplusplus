function r = convert_cells(mesh, convert_mode, varargin)
%CONVERT_CELLS  Linearize, simplexify or elevate the cells of a mesh.
%   R = mio.convert_cells(MESH, CONVERT_MODE, ...) with CONVERT_MODE
%   'linearize', 'simplexify', 'elevate' or 'reduce'. Options
%   'record_parent_ids' (attach 'convert:parent_cell'; false), 'order' (the
%   target order of 'elevate' or 'reduce'; 0 for 2 and 1) and 'full'
%   (full-Lagrange family when elevating; false). R is a struct: mesh, point_map
%   (n x 1, 1-based; 0 marks a pruned point) and cell_maps (one 1-based
%   n x 1 map per input block).

    o = mio.internal.opts(varargin, {'record_parent_ids', false, 'order', 0, 'full', false}, ...
        'mio.convert_cells');
    r = meshioplusplus_mex('mio_convert_cells', mesh, ...
        mio.internal.str(convert_mode, 'convert_mode'), logical(o.record_parent_ids), ...
        double(o.order), logical(o.full));
end
