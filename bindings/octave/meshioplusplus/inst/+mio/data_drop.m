function out = data_drop(mesh, location, names, varargin)
%DATA_DROP  A copy of MESH without the named data arrays at LOCATION.
%   OUT = mio.data_drop(MESH, LOCATION, NAMES) removes the arrays NAMES (a
%   string or cellstr; [] or {} drops nothing) at LOCATION ('point', 'cell'
%   or 'field'). Option: 'ignore_missing' (false) skips names that do not
%   exist instead of raising.
%
%   The data operations act on point, cell and field data only and never
%   touch geometry: points, connectivity and block order come through
%   bit-identical. The C++ data_manage is not exposed across the C ABI;
%   mio.data_drop, mio.data_keep and mio.data_rename compose to its effect.

    o = mio.internal.opts(varargin, {'ignore_missing', false}, 'mio.data_drop');
    out = meshioplusplus_mex('mio_data_drop', mesh, ...
        mio.internal.choice(location, {'point', 'cell', 'field'}, 'location'), ...
        mio.internal.cellstr_arg(names, 'names'), logical(o.ignore_missing));
end
