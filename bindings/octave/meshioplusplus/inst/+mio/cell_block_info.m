function out = cell_block_info(mesh, block)
%CELL_BLOCK_INFO  Shape of cell block BLOCK (1-based): a struct with num_cells,
%   nodes_per_cell, is_ragged, is_polyhedron, num_faces and num_nodes.

    out = meshioplusplus_mex('mio_cell_block_info', mesh, mio.internal.block(block));
end
