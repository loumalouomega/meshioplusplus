function m = interpolate(source, target, varargin)
%INTERPOLATE  Transfer point data from SOURCE onto the points of TARGET.
%   M = mio.interpolate(SOURCE, TARGET, ...) returns a new mesh: TARGET with
%   the transferred arrays. Options:
%     'method'        'nearest' or 'barycentric' ('nearest')
%     'arrays'        cellstr of source point-data names; [] transfers every
%                     source point-data array ([])
%     'extrapolate'   give a target point outside the source domain the
%                     nearest source value instead of 'default_value' (false)
%     'default_value' fill value outside the source domain (0)
%     'on_conflict'   'error', 'overwrite' or 'suffix' ('error')

    o = mio.internal.opts(varargin, {'method', 'nearest', 'arrays', [], ...
        'extrapolate', false, 'default_value', 0, 'on_conflict', 'error'}, 'mio.interpolate');
    m = meshioplusplus_mex('mio_interpolate', source, target, ...
        mio.internal.str(o.method, 'method'), mio.internal.cellstr_arg(o.arrays, 'arrays'), ...
        logical(o.extrapolate), double(o.default_value), ...
        mio.internal.str(o.on_conflict, 'on_conflict'));
end
