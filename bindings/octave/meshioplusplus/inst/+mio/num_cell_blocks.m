function out = num_cell_blocks(mesh)
%NUM_CELL_BLOCKS  Number of cell blocks.

    out = meshioplusplus_mex('mio_num_cell_blocks', mesh);
end
