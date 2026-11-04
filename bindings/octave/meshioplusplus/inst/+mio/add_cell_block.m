function mesh = add_cell_block(mesh, cell_type, conn)
%ADD_CELL_BLOCK  Append a cell block: CONN is num_cells x nodes_per_cell, 1-based.

    meshioplusplus_mex('mio_add_cell_block', mesh, mio.internal.str(cell_type, 'cell_type'), conn);
end
