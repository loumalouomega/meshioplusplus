function r = decimate(mesh, varargin)
%DECIMATE  Quadric-error decimation of a triangle surface.
%   R = mio.decimate(MESH, ...). Options:
%     'ratio'             fraction of faces to KEEP, in (0, 1]; negative = unset (-1)
%     'target_faces'      absolute face count to stop at; negative = unset (-1)
%     'max_error'         quadric-error ceiling; negative = unset (-1)
%     'placement'         'optimal', 'midpoint' or 'endpoint' ('optimal')
%     'preserve_boundary' pin boundary vertices (true)
%     'preserve_features' pin vertices on feature edges (true)
%     'feature_angle'     feature angle in degrees (30)
%     'frozen'            1-based point ids to pin outright; an out-of-range
%                         id raises by name ([])
%     'max_deviation'     the largest two-sided distance from the input on a
%                         'coverage' share of samples: the lightest mesh
%                         within it (roadmap 7.8.5); negative = unset (-1)
%     'coverage'          that share (0.999)
%     'max_faces'         at most this many faces, a cap on any criterion
%                         (alone, the target); negative = unset (-1)
%   R is a struct: mesh, point_map (n x 1, 1-based; 0 marks a removed point),
%   cell_maps, faces_removed, points_removed, collapses_rejected,
%   max_error_applied, deviation, deviation_met and rounds.

    o = mio.internal.opts(varargin, {'ratio', -1, 'target_faces', -1, 'max_error', -1, ...
        'placement', 'optimal', 'preserve_boundary', true, 'preserve_features', true, ...
        'feature_angle', 30, 'frozen', [], 'max_deviation', -1, 'coverage', 0.999, ...
        'max_faces', -1}, 'mio.decimate');
    r = meshioplusplus_mex('mio_decimate', mesh, double(o.ratio), double(o.target_faces), ...
        double(o.max_error), mio.internal.str(o.placement, 'placement'), ...
        logical(o.preserve_boundary), logical(o.preserve_features), ...
        double(o.feature_angle), double(o.frozen), double(o.max_deviation), ...
        double(o.coverage), double(o.max_faces));
end
