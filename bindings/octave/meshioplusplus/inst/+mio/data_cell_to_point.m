function out = data_cell_to_point(mesh, varargin)
%DATA_CELL_TO_POINT  Average cell data onto the points.
%   OUT = mio.data_cell_to_point(MESH) averages every cell data array onto
%   the points of a new mesh. Options:
%     'names'   arrays to convert ([] or {} for every cell data array)
%     'weight'  'uniform' (default) or 'measure' (weight by cell measure)
%     'suffix'  write to name + suffix instead of the same name ([])
%   Geometry is untouched (see mio.data_drop).

    o = mio.internal.opts(varargin, {'names', [], 'weight', 'uniform', 'suffix', []}, ...
        'mio.data_cell_to_point');
    out = meshioplusplus_mex('mio_data_cell_to_point', mesh, ...
        mio.internal.cellstr_arg(o.names, 'names'), ...
        mio.internal.choice(o.weight, {'uniform', 'measure'}, 'weight'), ...
        mio.internal.str(o.suffix, 'suffix'));
end
