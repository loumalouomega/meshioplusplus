function mesh = append_cell_data(mesh, name, data)
%APPEND_CELL_DATA  Append the array of cell data NAME for the next cell block (one row
%   per cell); call it once per block, in block order.

    meshioplusplus_mex('mio_append_cell_data', mesh, mio.internal.str(name, 'name'), data);
end
