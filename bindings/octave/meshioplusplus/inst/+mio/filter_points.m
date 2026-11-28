function r = filter_points(mesh, varargin)
%FILTER_POINTS  Remove outliers from, and thin, a point cloud.
%   R = mio.filter_points(MESH, ...) applies up to four steps in this order,
%   each seeing the points the one before kept: statistical outlier removal
%   ('statistical_k' neighbours, 'std_ratio' standard deviations above the mean;
%   Rusu et al. 2008), radius outlier removal ('radius', 'min_neighbors'), voxel
%   thinning ('voxel_size': one real point, the one nearest the mean, per cube)
%   and Poisson-disk thinning ('min_distance': kept points are at least that far
%   apart, greedy in index order). A step is off at 0 (the defaults, except
%   'std_ratio' 1 and 'min_neighbors' 4) and at least one must be on.
%   Returns a struct: mesh (the point cloud, per-point arrays gathered, Point
%   regions following, cells dropped), point_map (n x 1, each output point's
%   1-based input point), num_removed (1 x 4: statistical, radius, voxel,
%   poisson) and warnings (a cell array of char).

    o = mio.internal.opts(varargin, {'statistical_k', 0, 'std_ratio', 1, 'radius', 0, ...
        'min_neighbors', 4, 'voxel_size', 0, 'min_distance', 0, 'voxel_mode', 'nearest'}, ...
        'mio.filter_points');
    r = meshioplusplus_mex('mio_filter_points', mesh, double(o.statistical_k), ...
        double(o.std_ratio), double(o.radius), double(o.min_neighbors), ...
        double(o.voxel_size), double(o.min_distance), ...
        mio.internal.str(o.voxel_mode, 'voxel_mode'));
    if isempty(r.warnings)
        r.warnings = {};
    else
        r.warnings = strsplit(r.warnings, char(10));
    end
end
