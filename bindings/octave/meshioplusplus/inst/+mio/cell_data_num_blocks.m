function n = cell_data_num_blocks(mesh, name)
%CELL_DATA_NUM_BLOCKS  Number of per-block arrays of cell data NAME.

    n = meshioplusplus_mex('mio_cell_data_num_blocks', mesh, mio.internal.str(name, 'name'));
end
