function out = data_rename(mesh, location, from, to)
%DATA_RENAME  A copy of MESH with data array FROM at LOCATION renamed TO.
%   OUT = mio.data_rename(MESH, LOCATION, FROM, TO) keeps the values, class
%   and shape. LOCATION is 'point', 'cell' or 'field'. A missing FROM is an
%   error. Geometry is untouched (see mio.data_drop).

    out = meshioplusplus_mex('mio_data_rename', mesh, ...
        mio.internal.choice(location, {'point', 'cell', 'field'}, 'location'), ...
        mio.internal.str(from, 'from'), mio.internal.str(to, 'to'));
end
