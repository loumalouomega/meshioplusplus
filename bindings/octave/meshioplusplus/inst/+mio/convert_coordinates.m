function r = convert_coordinates(mesh, varargin)
%CONVERT_COORDINATES  Cylindrical or spherical components about an axis.
%   The points stay put; their coordinates are recorded as point data
%   ('coords:cylindrical' or 'coords:spherical'), and each named vector or
%   tensor array is rewritten in the local basis. Options: 'system'
%   ('cylindrical' or 'spherical'), 'axis' ([0 0 1]), 'origin' ([0 0 0]),
%   'reference' ([1 0 0]), 'arrays' ({}), 'inverse' (false) and
%   'record_coordinates' (true). Returns a struct: mesh, num_arrays and
%   num_singular.

    o = mio.internal.opts(varargin, {'system', 'cylindrical', 'axis', [0 0 1], ...
        'origin', [0 0 0], 'reference', [1 0 0], 'arrays', {}, 'inverse', false, ...
        'record_coordinates', true}, 'mio.convert_coordinates');
    r = meshioplusplus_mex('mio_convert_coordinates', mesh, ...
        mio.internal.str(o.system, 'system'), double(o.axis), double(o.origin), ...
        double(o.reference), mio.internal.cellstr_arg(o.arrays, 'arrays'), ...
        logical(o.inverse), logical(o.record_coordinates));
end
