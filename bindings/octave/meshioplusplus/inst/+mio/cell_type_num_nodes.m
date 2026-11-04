function n = cell_type_num_nodes(name)
%CELL_TYPE_NUM_NODES  Nodes per cell of cell type NAME (-1 for variable-size types).

    n = meshioplusplus_mex('mio_cell_type_num_nodes', mio.internal.str(name, 'name'));
end
