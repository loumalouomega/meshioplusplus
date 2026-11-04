function r = optimize_volume(mesh, varargin)
%OPTIMIZE_VOLUME  ODT-remesh a tetrahedral mesh (relocation and flips).
%   R = mio.optimize_volume(MESH, ...) raises the worst element quality of a
%   tet-only mesh by relocating vertices and changing connectivity. The point
%   set is invariant, so point data and point regions carry; cell data and
%   cell/side regions are dropped. Options: 'max_iterations' (10),
%   'relocate' (true), 'flip' (true), 'preserve_boundary' (true),
%   'min_improvement' (1e-6), 'objective' ('worst' or 'energy'; 'worst') and
%   'max_ring_size' (7, 'worst' only). Returns a struct: mesh, num_flips,
%   num_23_flips, num_32_flips, num_vertices_moved, num_tets,
%   min_quality_before, min_quality_after, num_edge_removals,
%   num_multi_face_removals, num_worst_smoothed, min_dihedral_before and
%   min_dihedral_after.

    o = mio.internal.opts(varargin, {'max_iterations', 10, 'relocate', true, 'flip', true, ...
        'preserve_boundary', true, 'min_improvement', 1e-6, 'objective', 'worst', ...
        'max_ring_size', 7}, 'mio.optimize_volume');
    r = meshioplusplus_mex('mio_optimize_volume', mesh, double(o.max_iterations), ...
        logical(o.relocate), logical(o.flip), logical(o.preserve_boundary), ...
        double(o.min_improvement), mio.internal.str(o.objective, 'objective'), ...
        double(o.max_ring_size));
end
