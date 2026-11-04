function t = cell_block_types(mesh)
%CELL_BLOCK_TYPES  The cell type of every block (cellstr, n x 1).

    nb = mio.num_cell_blocks(mesh);
    t = cell(nb, 1);
    for b = 1:nb
        t{b} = mio.cell_block_type(mesh, b);
    end
end
