function r = find_symmetry_plane(mesh, varargin)
%FIND_SYMMETRY_PLANE  The plane MESH is most nearly symmetric about.
%   Each principal plane through the centroid (and an optional seed) is
%   mirrored, registered back by ICP and read off the reflection-motion
%   product (roadmap 7.3.6). Options: 'seed' ({origin, normal}, {}),
%   'tolerance' (0: 0.005 of the diagonal), 'min_share' (0.9) and
%   'max_samples' (4000). Returns a struct: origin, normal, rms, share,
%   tolerance and chosen_candidate.

    o = mio.internal.opts(varargin, {'seed', {}, 'tolerance', 0, 'min_share', 0.9, ...
        'max_samples', 4000}, 'mio.find_symmetry_plane');
    so = [];
    sn = [];
    if ~isempty(o.seed)
        so = double(o.seed{1});
        sn = double(o.seed{2});
    end
    r = meshioplusplus_mex('mio_find_symmetry_plane', mesh, so, sn, double(o.tolerance), ...
        double(o.min_share), double(o.max_samples));
end
