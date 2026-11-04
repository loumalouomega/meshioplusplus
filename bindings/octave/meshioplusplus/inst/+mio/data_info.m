function info = data_info(mesh)
%DATA_INFO  Per-array summaries of every data array.
%   INFO = mio.data_info(MESH) returns an n x 1 struct array, one entry per
%   point, cell and field data array, with fields name, location ('point',
%   'cell' or 'field'), dtype, num_blocks, num_entries, num_components,
%   num_values, min, max, mean, num_nan, num_inf, num_finite,
%   inconsistent_blocks (logical) and components: a (num_components x 3)
%   matrix whose columns are the per-component min, max and mean.

    info = meshioplusplus_mex('mio_data_info', mesh);
end
