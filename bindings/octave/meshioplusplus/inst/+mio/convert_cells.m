function r = convert_cells(mesh, convert_mode, varargin)
%CONVERT_CELLS  Linearize, simplexify or elevate the cells of a mesh.
%   R = mio.convert_cells(MESH, CONVERT_MODE, ...) with CONVERT_MODE
%   'linearize', 'simplexify' or 'elevate'. Option 'record_parent_ids'
%   (attach 'convert:parent_cell'; false). R is a struct: mesh, point_map
%   (n x 1, 1-based; 0 marks a pruned point) and cell_maps (one 1-based
%   n x 1 map per input block).

    o = mio.internal.opts(varargin, {'record_parent_ids', false}, 'mio.convert_cells');
    r = meshioplusplus_mex('mio_convert_cells', mesh, ...
        mio.internal.str(convert_mode, 'convert_mode'), logical(o.record_parent_ids));
end
