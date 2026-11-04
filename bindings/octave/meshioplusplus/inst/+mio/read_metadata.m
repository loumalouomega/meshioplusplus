function s = read_metadata(path, varargin)
%READ_METADATA  Summarise a mesh file without reading its arrays.
%   Returns a struct: num_points, point_dim, num_cells, num_cell_blocks,
%   cell_block_types, cell_block_num_cells, cell_block_nodes_per_cell,
%   cell_block_is_ragged, point/cell/field_data_names, bbox (2 x 3: min row,
%   max row; [] when unknown), fell_back, time_values and regions (struct
%   array without entries). Option: 'format'.

    o = mio.internal.opts(varargin, {'format', []}, 'mio.read_metadata');
    s = meshioplusplus_mex('mio_read_metadata', mio.internal.str(path, 'path'), ...
        mio.internal.str(o.format, 'format'));
end
