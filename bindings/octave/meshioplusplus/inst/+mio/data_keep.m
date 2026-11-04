function out = data_keep(mesh, location, names, varargin)
%DATA_KEEP  A copy of MESH keeping only the named data arrays at LOCATION.
%   OUT = mio.data_keep(MESH, LOCATION, NAMES) keeps the arrays NAMES (a
%   string or cellstr) at LOCATION ('point', 'cell' or 'field') and drops the
%   rest there. Option: 'ignore_missing' (false) skips names that do not
%   exist. Geometry is untouched (see mio.data_drop).

    o = mio.internal.opts(varargin, {'ignore_missing', false}, 'mio.data_keep');
    out = meshioplusplus_mex('mio_data_keep', mesh, ...
        mio.internal.choice(location, {'point', 'cell', 'field'}, 'location'), ...
        mio.internal.cellstr_arg(names, 'names'), logical(o.ignore_missing));
end
