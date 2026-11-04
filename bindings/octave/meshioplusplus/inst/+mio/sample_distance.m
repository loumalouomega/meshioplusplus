function d = sample_distance(mesh, points, varargin)
%SAMPLE_DISTANCE  Signed distance from a surface at arbitrary query points.
%   D = mio.sample_distance(MESH, POINTS) evaluates the distance from the
%   surface MESH at POINTS, an (n x 3) matrix with one query point per row,
%   and returns an n x 1 column (negative inside). Options:
%     'sign'              'pseudonormal' (default), 'winding-number' or
%                         'unsigned'
%     'band'              clamp distances beyond this; 0 is the full field (0)
%     'watertight_check'  'off', 'warn' (default) or 'error'

    o = mio.internal.opts(varargin, {'sign', 'pseudonormal', 'band', 0, ...
        'watertight_check', 'warn'}, 'mio.sample_distance');
    sign = mio.internal.od_code(o.sign, {'unsigned', 'pseudonormal', 'winding-number'}, 'sign');
    check = mio.internal.od_code(o.watertight_check, {'off', 'warn', 'error'}, 'watertight check');
    d = meshioplusplus_mex('mio_sample_distance', mesh, double(points), sign, ...
        double(o.band), check);
end
