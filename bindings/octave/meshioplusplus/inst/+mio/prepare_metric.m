function r = prepare_metric(mesh, arrays, varargin)
%PREPARE_METRIC  Intersect, normalise, clamp and grade metric fields.
%   ARRAYS names the metric point data (6 components in 3-D, 3 in 2-D).
%   Options: 'output' ('metric'), 'target_vertices' (0: skip), 'norm' (2; 0
%   is infinity), 'hmin' and 'hmax' (0: that side open) and 'gradation' (0:
%   skip). Returns a struct: mesh, complexity_before, complexity_after and
%   num_graded. See doc/metric.md.

    o = mio.internal.opts(varargin, {'output', 'metric', 'target_vertices', 0, 'norm', 2, ...
        'hmin', 0, 'hmax', 0, 'gradation', 0}, 'mio.prepare_metric');
    r = meshioplusplus_mex('mio_prepare_metric', mesh, mio.internal.cellstr_arg(arrays, 'arrays'), ...
        mio.internal.str(o.output, 'output'), double(o.target_vertices), double(o.norm), ...
        double(o.hmin), double(o.hmax), double(o.gradation));
end
