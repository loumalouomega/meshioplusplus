function r = adapt(mesh, varargin)
%ADAPT  Remesh to a size or a metric with MMG.
%   R = mio.adapt(MESH, ...) remeshes through mmg3d for tetrahedra, mmg2d for
%   a planar triangle mesh, mmgs for another triangle surface (roadmap 3.2.1,
%   doc/mmg.md). Give at most one of 'hsiz', 'size' and 'metric'. Cell and
%   side regions come back under their names. Fails naming
%   MESHIOPLUSPLUS_WITH_MMG without the backend (mio.has_mmg). Options:
%     'hsiz'            a constant target edge size, negative for none (-1)
%     'size'            a one-component point-data name: the size per point ('')
%     'metric'          a six-component (three in 2-D) point-data name: the
%                       metric in the lower-triangle order of a Medit .sol ('')
%     'hmin', 'hmax', 'hausd', 'hgrad'
%                       MMG's sizes, Hausdorff distance and gradation,
%                       negative for MMG's default (-1 each)
%     'angle_detection' the ridge angle in degrees; 0 turns detection off (45)
%     'frozen'          1-based points MMG must keep ([])
%     'required'        cellstr of regions whose entities MMG must keep ({})
%   R is a struct: mesh, status ('success' or 'low_failure'), library,
%   num_points_in, num_points_out, num_elements_in, num_elements_out, quality
%   (3 x 1: min, mean, max of 'quality:aspect_ratio') and num_warnings.

    o = mio.internal.opts(varargin, {'hsiz', -1, 'size', '', 'metric', '', 'hmin', -1, ...
        'hmax', -1, 'hausd', -1, 'hgrad', -1, 'angle_detection', 45, 'frozen', [], ...
        'required', {}}, 'mio.adapt');
    r = meshioplusplus_mex('mio_adapt', mesh, double(o.hsiz), mio.internal.str(o.size, 'size'), ...
        mio.internal.str(o.metric, 'metric'), double(o.hmin), double(o.hmax), ...
        double(o.hausd), double(o.hgrad), double(o.angle_detection), double(o.frozen), ...
        mio.internal.cellstr_arg(o.required, 'required'));
end
