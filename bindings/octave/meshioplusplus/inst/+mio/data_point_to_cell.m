function out = data_point_to_cell(mesh, varargin)
%DATA_POINT_TO_CELL  Average point data onto the cells.
%   OUT = mio.data_point_to_cell(MESH) averages every point data array onto
%   the cells of a new mesh. Options:
%     'names'   arrays to convert ([] or {} for every point data array)
%     'suffix'  write to name + suffix instead of the same name ([])
%   Geometry is untouched (see mio.data_drop).

    o = mio.internal.opts(varargin, {'names', [], 'suffix', []}, 'mio.data_point_to_cell');
    out = meshioplusplus_mex('mio_data_point_to_cell', mesh, ...
        mio.internal.cellstr_arg(o.names, 'names'), mio.internal.str(o.suffix, 'suffix'));
end
