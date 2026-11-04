function mesh = add_polygon_block(mesh, cell_type, rows)
%ADD_POLYGON_BLOCK  Append a ragged polygon block: ROWS is a cell array of 1-based node
%   vectors, one per cell. Returns MESH.

    meshioplusplus_mex('mio_add_polygon_block', mesh, mio.internal.str(cell_type, 'cell_type'), rows);
end
