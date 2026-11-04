function n = num_cells(mesh)
%NUM_CELLS  Total number of cells over every block.

    n = 0;
    for b = 1:mio.num_cell_blocks(mesh)
        info = mio.cell_block_info(mesh, b);
        n = n + info.num_cells;
    end
end
