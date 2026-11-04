function r = repair(mesh, varargin)
%REPAIR  Repair a surface's orientation, holes and pinched vertices.
%   Options: 'fix_orientation' (true), 'orient_outward' (true), 'fill_holes'
%   (true), 'split_non_manifold' (true), 'max_hole_edges' (10; <= 0 no
%   limit), 'weld_tolerance' (0: no weld) and 'record_provenance' (false).
%   Returns a struct: mesh, quality_before, quality_after (structs of
%   boundary_edges, non_manifold_edges, inconsistent_pairs,
%   degenerate_triangles, watertight), num_flipped, num_components,
%   largest_component, num_oriented_outward, num_unorientable,
%   num_vertices_split, num_holes_detected, num_holes_filled,
%   num_holes_skipped, num_faces_added, num_points_added and points_welded.

    o = mio.internal.opts(varargin, {'fix_orientation', true, 'orient_outward', true, ...
        'fill_holes', true, 'split_non_manifold', true, 'max_hole_edges', 10, ...
        'weld_tolerance', 0, 'record_provenance', false}, 'mio.repair');
    r = meshioplusplus_mex('mio_repair', mesh, logical(o.fix_orientation), ...
        logical(o.orient_outward), logical(o.fill_holes), logical(o.split_non_manifold), ...
        double(o.max_hole_edges), double(o.weld_tolerance), logical(o.record_provenance));
end
