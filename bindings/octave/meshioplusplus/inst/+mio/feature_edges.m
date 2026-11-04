function r = feature_edges(mesh, varargin)
%FEATURE_EDGES  Sharp, open, non-manifold and inconsistent edges as a line mesh.
%   Options: 'feature_angle' (30 degrees, in [0, 180]), 'feature',
%   'boundary', 'non_manifold', 'inconsistent' (all true: report that
%   category) and 'region' (''). Returns a struct: mesh (with cell data
%   'feature:kind' -- 1 feature, 2 boundary, 3 non-manifold, 4 inconsistent
%   -- and 'feature:angle'), num_feature, num_boundary, num_non_manifold and
%   num_inconsistent.

    o = mio.internal.opts(varargin, {'feature_angle', 30, 'feature', true, 'boundary', true, ...
        'non_manifold', true, 'inconsistent', true, 'region', ''}, 'mio.feature_edges');
    r = meshioplusplus_mex('mio_feature_edges', mesh, double(o.feature_angle), ...
        logical(o.feature), logical(o.boundary), logical(o.non_manifold), ...
        logical(o.inconsistent), mio.internal.str(o.region, 'region'));
end
