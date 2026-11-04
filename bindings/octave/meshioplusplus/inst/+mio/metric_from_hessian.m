function m = metric_from_hessian(mesh, hessian, varargin)
%METRIC_FROM_HESSIAN  An anisotropic metric from a Hessian.
%   M = mio.metric_from_hessian(MESH, HESSIAN, ...) computes the metric of
%   the row-major point-data Hessian HESSIAN ((n x 9), or (n x 4) in 2-D)
%   for mio.adapt: M = c/error |H|, its sizes clamped to [hmin, hmax], as
%   point data 'output'. Needs no MMG. Options:
%     'error'  the target interpolation error (1e-2)
%     'hmin'   the smallest edge size the metric asks for (1e-3)
%     'hmax'   the largest edge size the metric asks for (1)
%     'output' the result's point-data name ('metric')
%   M is a new mesh with the metric.

    o = mio.internal.opts(varargin, {'error', 1e-2, 'hmin', 1e-3, 'hmax', 1, ...
        'output', 'metric'}, 'mio.metric_from_hessian');
    m = meshioplusplus_mex('mio_metric_from_hessian', mesh, ...
        mio.internal.str(hessian, 'hessian'), double(o.error), double(o.hmin), ...
        double(o.hmax), mio.internal.str(o.output, 'output'));
end
