function r = mesh_cad_surface(info, target_size, varargin)
%MESH_CAD_SURFACE  Triangulate a CAD model's faces.
%   R = mio.mesh_cad_surface(INFO, TARGET_SIZE, ...) builds a conforming
%   triangle mesh of the faces of the B-rep the mio.FormatInfo INFO holds
%   (Kratos CAD JSON, STEP or IGES, from mio.read_with_info) with edges of
%   about TARGET_SIZE (> 0), with no CAD kernel (v17.13.0). Every edge is
%   sampled once and shared by the faces that use it, so a closed solid comes
%   out watertight; each face is a constrained Delaunay triangulation of its
%   parameters lifted onto its surface. A face that cannot be meshed is left
%   out and reported in face_status / face_reason, never meshed wrong;
%   'strict' makes that an error. See doc/mesh_cad.md.
%
%   Options (positional in this order, or name/value; R defaults):
%     'max_chordal_deviation'  sagitta cap of boundary segments; 0 off (0)
%     'min_points_per_edge'    at least this many points inside an edge (0)
%     'tolerance'              geometric tolerance; 0 is 1e-7 of the size (0)
%     'record_face_ids'        'cad:face' cell data and face:<k> regions (true)
%     'edges'                  a line block of the edges, 'cad:edge' (true)
%     'match_unreferenced_trims'  pair trims no edge names (true)
%     'strict'                 fail when any face cannot be meshed (false)
%     'max_points'             the most points in the result (50000000)
%     'min_angle'              refine to this smallest angle, degrees < 60;
%                              0 keeps the unrefined mesh (0)
%     'smooth_iterations'      guarded smoothing passes (0)
%     'smooth_method'          'laplacian' or 'odt'
%     'order'                  1 linear; 2 or 3 curved triangle6/triangle10
%                              (and line3/line4) on the geometry (1)
%     'grading'                largest size growth; 0 is 1.5 (0)
%     'curvature_deviation'    chord-to-surface bound; 0 off (0)
%     'cells_across_gaps'      elements across the local feature size (0)
%     'pole_rings'             rings seeded around each pole, <= 8 (0)
%     'face_sizes'             (n x 2) [1-based face index, size] rows, or []
%     'edge_sizes'             the same per edge, or []
%     'background'             a mesh whose point data 'background_array' is
%                              the size, or []
%     'background_array'       ('size')
%     'metric'                 'mean' or 'varying'
%     'virtual_edges'          sample chains of short smooth edges as one
%                              (false)
%     'reparam_condition_number'  flatten a face whose metric's condition
%                              number passes this; 0 off, else > 1 (0)
%     'model_tolerances'       accept the gaps the model declares (false)
%     'virtual_faces'          mesh a narrow face with its neighbour (false)
%     'lop_corners'            lop acute corners; only with refinement
%                              (false)
%     'metric_background'      a mesh with a 9-component Hessian (or the
%                              1-component field) 'metric_array', or []
%     'metric_array'           ('')
%     'metric_error'           target interpolation error (0)
%     'untangle'               with order 2 or 3, untangle folded curved
%                              elements (false)
%     'quads'                  'none', 'dominant' or 'all'
%     'quad_pairs'             'matching' or 'greedy'
%
%   R is a struct: mesh; the counts num_points, num_triangles,
%   num_edge_segments, num_faces, num_failed, num_edges, num_open_edges,
%   num_non_manifold_edges, num_misoriented_edges, num_degenerate_dropped;
%   tolerance, min_angle_deg, max_angle_deg, min_edge, max_edge,
%   max_edge_gap; face_status and face_reason (n x 1 cellstr, one per face);
%   warnings (cellstr); num_refinement_points, num_faces_below_bound,
%   num_smoothed, min_input_angle_deg; angle_histogram (18 x 1, ten-degree
%   bins); corners; num_curved_straightened; face_size_limits (cell, per face
%   the 1-based edges finer than its size allows); virtual_edges (cell, each
%   chain's 1-based edges); reparametrised, merged_faces, tolerance_uses,
%   lopped_corners; num_untangled, num_quads, num_leftover_triangles,
%   min_quad_scaled_jacobian and min_metric_angle.
%
%   The tables (R data frames) corners (face, point, u, v, x, y, z,
%   angle_deg, min_forced_deg, reach, pole), reparametrised (face,
%   condition_number, condition_number_disc, method: 1 conformal, 2
%   mean-value, 3 annulus, 0 declined, reason), merged_faces (sliver, host,
%   width, merged, reason), tolerance_uses (edge, declared, gap, where) and
%   lopped_corners (face, point, angle_deg, lfs0, lop_radius, turning_deg,
%   stated_bound_deg, num_forced, min_forced_deg, within_two_radii,
%   guaranteed, reason) are structs of n x 1 columns: R.corners.angle_deg is
%   a column and numel(R.corners.face) the row count. Faces, points and edges
%   in them are 1-based.

    o = mio.internal.opts(varargin, {'max_chordal_deviation', 0, 'min_points_per_edge', 0, ...
        'tolerance', 0, 'record_face_ids', true, 'edges', true, ...
        'match_unreferenced_trims', true, 'strict', false, 'max_points', 50000000, ...
        'min_angle', 0, 'smooth_iterations', 0, 'smooth_method', 'laplacian', 'order', 1, ...
        'grading', 0, 'curvature_deviation', 0, 'cells_across_gaps', 0, 'pole_rings', 0, ...
        'face_sizes', [], 'edge_sizes', [], 'background', [], 'background_array', 'size', ...
        'metric', 'mean', 'virtual_edges', false, 'reparam_condition_number', 0, ...
        'model_tolerances', false, 'virtual_faces', false, 'lop_corners', false, ...
        'metric_background', [], 'metric_array', '', 'metric_error', 0, 'untangle', false, ...
        'quads', 'none', 'quad_pairs', 'matching'}, 'mio.mesh_cad_surface');
    smooth_method = mio.internal.choice(o.smooth_method, {'laplacian', 'odt'}, 'smooth_method');
    metric = mio.internal.choice(o.metric, {'mean', 'varying'}, 'metric');
    quads = mio.internal.choice(o.quads, {'none', 'dominant', 'all'}, 'quads');
    quad_pairs = mio.internal.choice(o.quad_pairs, {'matching', 'greedy'}, 'quad_pairs');
    [fidx, fsz] = mio.internal.fi_size_pairs(o.face_sizes, 'face_sizes');
    [eidx, esz] = mio.internal.fi_size_pairs(o.edge_sizes, 'edge_sizes');
    r = meshioplusplus_mex('mio_mesh_cad_surface', info, double(target_size), ...
        double(o.max_chordal_deviation), double(o.min_points_per_edge), double(o.tolerance), ...
        logical(o.record_face_ids), logical(o.edges), logical(o.match_unreferenced_trims), ...
        logical(o.strict), double(o.max_points), double(o.min_angle), ...
        double(o.smooth_iterations), smooth_method, double(o.order), double(o.grading), ...
        double(o.curvature_deviation), double(o.cells_across_gaps), double(o.pole_rings), ...
        fidx, fsz, eidx, esz, o.background, mio.internal.str(o.background_array, 'background_array'), ...
        metric, logical(o.virtual_edges), double(o.reparam_condition_number), ...
        logical(o.model_tolerances), logical(o.virtual_faces), logical(o.lop_corners), ...
        o.metric_background, mio.internal.str(o.metric_array, 'metric_array'), ...
        double(o.metric_error), logical(o.untangle), quads, quad_pairs);
    r.warnings = mio.internal.fi_split_lines(r.warnings);
end
