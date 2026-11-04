function out = cell_data_names(mesh)
%CELL_DATA_NAMES  Cell-data names, sorted (cellstr).

    out = meshioplusplus_mex('mio_cell_data_names', mesh);
end
