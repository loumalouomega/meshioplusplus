function r = gradient(mesh, array, varargin)
%GRADIENT  Gradient, divergence or curl of the point data ARRAY.
%   R = mio.gradient(MESH, ARRAY, ...) with options 'op' ('gradient',
%   'divergence' or 'curl'), 'method' ('green-gauss' or 'least-squares'),
%   'location' ('cell' or 'point'), 'output' ('' uses '<array>:<op>'),
%   'component' (-1: every component) and 'overwrite' (false). Returns a
%   struct: mesh, num_skipped and num_fallback.

    o = mio.internal.opts(varargin, {'op', 'gradient', 'method', 'green-gauss', ...
        'location', 'cell', 'output', '', 'component', -1, 'overwrite', false}, 'mio.gradient');
    r = meshioplusplus_mex('mio_gradient', mesh, mio.internal.str(array, 'array'), ...
        mio.internal.str(o.op, 'op'), mio.internal.str(o.method, 'method'), ...
        mio.internal.str(o.location, 'location'), mio.internal.str(o.output, 'output'), ...
        double(o.component), logical(o.overwrite));
end
