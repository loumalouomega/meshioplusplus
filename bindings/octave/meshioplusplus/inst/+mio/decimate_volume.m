function r = decimate_volume(mesh, varargin)
%DECIMATE_VOLUME  Quadric-error decimation of a tetrahedral mesh.
%   R = mio.decimate_volume(MESH, ...) collapses tet edges greedily by
%   boundary quadric error; exactly one of 'ratio', 'target_cells' and
%   'max_error' must be set. MESH must be tet-only (simplexify it with
%   mio.convert_cells first). Options:
%     'ratio'             fraction of tets to KEEP, in (0, 1]; negative = unset (-1)
%     'target_cells'      absolute tet count to stop at; negative = unset (-1)
%     'max_error'         boundary quadric-error ceiling; negative = unset (-1)
%     'placement'         'optimal', 'midpoint' or 'endpoint' ('optimal')
%     'preserve_boundary' pin every boundary vertex (false; mio.decimate's
%                         default is true)
%     'preserve_features' pin boundary vertices on feature edges (true)
%     'feature_angle'     feature angle in degrees (30)
%     'frozen'            1-based point ids to pin outright; an out-of-range
%                         id raises by name ([])
%   R is a struct: mesh, point_map (n x 1, 1-based; 0 marks a removed point),
%   cell_maps (one n x 1 column per input block; 0 marks a removed tet),
%   tets_removed, points_removed, collapses_rejected and max_error_applied.

    o = mio.internal.opts(varargin, {'ratio', -1, 'target_cells', -1, 'max_error', -1, ...
        'placement', 'optimal', 'preserve_boundary', false, 'preserve_features', true, ...
        'feature_angle', 30, 'frozen', []}, 'mio.decimate_volume');
    r = meshioplusplus_mex('mio_decimate_volume', mesh, double(o.ratio), ...
        double(o.target_cells), double(o.max_error), ...
        mio.internal.str(o.placement, 'placement'), logical(o.preserve_boundary), ...
        logical(o.preserve_features), double(o.feature_angle), double(o.frozen));
end
