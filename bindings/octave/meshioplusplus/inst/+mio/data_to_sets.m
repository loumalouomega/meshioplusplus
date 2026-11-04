function out = data_to_sets(mesh, location, key)
%DATA_TO_SETS  Convert one scalar integer data array into region-backed sets.
%   OUT = mio.data_to_sets(MESH, LOCATION, KEY) returns a new mesh with one
%   set per label of the array KEY at LOCATION ('point' or 'cell'); KEY
%   itself is removed. The inverse of mio.sets_to_data.

    out = meshioplusplus_mex('mio_data_to_sets', mesh, ...
        mio.internal.choice(location, {'point', 'cell', 'field'}, 'location'), ...
        mio.internal.str(key, 'key'));
end
