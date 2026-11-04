function out = cell_block_type(mesh, block)
%CELL_BLOCK_TYPE  The cell type name of block BLOCK (1-based).

    out = meshioplusplus_mex('mio_cell_block_type', mesh, mio.internal.block(block));
end
