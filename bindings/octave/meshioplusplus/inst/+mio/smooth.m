function r = smooth(mesh, varargin)
%SMOOTH  Move the free nodes of MESH by Taubin, Laplacian or ODT smoothing.
%   R = mio.smooth(MESH, ...) with options:
%     'method'            'taubin' (default), 'laplacian' or 'odt' (tet-only)
%     'iterations'        smoothing iterations (10)
%     'lambda'            relaxation factor in (0, 1); negative (-1) is the
%                         method's own default (0.5 laplacian, 0.33 taubin,
%                         0.9 odt)
%     'mu'                Taubin's un-shrinking factor, mu < -lambda < 0 (-0.34)
%     'fix_boundary'      pin every node on a boundary facet (true)
%     'preserve_features' pin boundary nodes at creases above
%                         'feature_angle' (true)
%     'feature_angle'     degrees (30)
%     'guard_inversion'   reject moves that flip a cell (true)
%     'frozen'            1-based point ids to pin outright ([] for none)
%   Returns a struct: mesh, num_nodes_moved, max_displacement and
%   num_skipped_inversion.

    o = mio.internal.opts(varargin, {'method', 'taubin', 'iterations', 10, 'lambda', -1, ...
        'mu', -0.34, 'fix_boundary', true, 'preserve_features', true, 'feature_angle', 30, ...
        'guard_inversion', true, 'frozen', []}, 'mio.smooth');
    r = meshioplusplus_mex('mio_smooth', mesh, mio.internal.str(o.method, 'method'), ...
        double(o.iterations), double(o.lambda), double(o.mu), logical(o.fix_boundary), ...
        logical(o.preserve_features), double(o.feature_angle), logical(o.guard_inversion), ...
        double(o.frozen));
end
