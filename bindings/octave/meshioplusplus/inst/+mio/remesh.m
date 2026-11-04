function r = remesh(mesh, num_clusters, varargin)
%REMESH  Remesh a surface by ACVD clustering at NUM_CLUSTERS vertices.
%   The output has new points and connectivity: point/cell data and regions
%   are dropped, field data is carried. Options: 'subdivide' (-1: automatic,
%   0 off), 'subsample_ratio' (10), 'max_subdivide' (4), 'max_iterations'
%   (100), 'max_repair_passes' (10), 'metric' ('isotropic', 'quadric' or
%   'anisotropic'), 'gradation' (0), 'preserve_boundary' (true) and
%   'max_anisotropy' (4, 'anisotropic' only). Returns a struct: mesh,
%   num_clusters, num_iterations, subdivide_applied, num_isolated_clusters
%   and num_non_manifold_vertices.

    o = mio.internal.opts(varargin, {'subdivide', -1, 'subsample_ratio', 10, ...
        'max_subdivide', 4, 'max_iterations', 100, 'max_repair_passes', 10, ...
        'metric', 'isotropic', 'gradation', 0, 'preserve_boundary', true, ...
        'max_anisotropy', 4}, 'mio.remesh');
    r = meshioplusplus_mex('mio_remesh', mesh, double(num_clusters), double(o.subdivide), ...
        double(o.subsample_ratio), double(o.max_subdivide), double(o.max_iterations), ...
        double(o.max_repair_passes), mio.internal.str(o.metric, 'metric'), ...
        double(o.gradation), logical(o.preserve_boundary), double(o.max_anisotropy));
end
