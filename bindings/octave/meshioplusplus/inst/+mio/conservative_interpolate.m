function m = conservative_interpolate(source, target, varargin)
%CONSERVATIVE_INTERPOLATE  Conservative (integral-preserving) data transfer.
%   M = mio.conservative_interpolate(SOURCE, TARGET, ...) returns a new mesh:
%   TARGET with the transferred arrays. Options:
%     'arrays'        cellstr of source array names; [] transfers every array ([])
%     'default_value' fill value outside the source domain (0)
%     'on_conflict'   'error', 'overwrite' or 'suffix' ('error')

    o = mio.internal.opts(varargin, {'arrays', [], 'default_value', 0, ...
        'on_conflict', 'error'}, 'mio.conservative_interpolate');
    m = meshioplusplus_mex('mio_conservative_interpolate', source, target, ...
        mio.internal.cellstr_arg(o.arrays, 'arrays'), double(o.default_value), ...
        mio.internal.str(o.on_conflict, 'on_conflict'));
end
