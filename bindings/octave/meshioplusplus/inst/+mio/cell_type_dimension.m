function n = cell_type_dimension(name)
%CELL_TYPE_DIMENSION  Topological dimension of cell type NAME.

    n = meshioplusplus_mex('mio_cell_type_dimension', mio.internal.str(name, 'name'));
end
