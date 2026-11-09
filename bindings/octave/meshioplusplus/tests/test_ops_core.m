function tests = test_ops_core()
%TEST_OPS_CORE  The operations of ops_core.c: surfaces, transforms, cleaning,
%   subsetting, derivatives, segmentation and CAD fitting, interfaces, repair
%   and the geometry kernels. Ported from the R suite, with smoke tests for
%   the functions it does not exercise.
    tests = {
        'geometry operations work', @t_geometry
        'transform moves points and leaves topology alone', @t_transform
        'clean, smooth and crop work', @t_clean_smooth_crop
        'slice and isosurface work', @t_slice_iso
        'gradient differentiates a point-data field', @t_gradient
        'hessian computes the second derivative of a scalar field', @t_hessian
        'estimate_error estimates a recovery-based indicator and marks', @t_estimate_error
        'crop_predicate keeps the cells a data comparison selects', @t_crop_predicate
        'remesh replaces a surface''s triangulation by ACVD clustering', @t_remesh
        'compute_normals gives smooth or split normals', @t_normals
        'compute_curvature totals the angle defect of a closed surface', @t_curvature
        'repair rewinds, fills and splits a surface', @t_repair
        'shrinkwrap projects points onto a target surface', @t_shrinkwrap
        'sobolev_deform filters a displacement field', @t_sobolev
        'feature edges, segmentation and primitive fits', @t_feature_segment
        'refine_patches refines a segmentation', @t_refine_patches
        'fit_cad builds a solid or falls back to native', @t_fit_cad
        'Hausdorff distance and the quality gate', @t_hausdorff_quality
        'geodesic distance, atlas, morph, untangle and metric', @t_geodesic_etc
        'eigenmodes, registration and point normals', @t_spectral_registration
        'reconstruct_surface closes a sampled sphere', @t_reconstruct
        'hulls, Delaunay and coordinate conversion', @t_hull_coords
        'interface search, contact projection and interface splitting', @t_interface
        'region adjacency between two cell regions', @t_region_adjacency
    };
end

% --- shared builders --------------------------------------------------------

function m = tri_grid(n)
% An n x n point grid in the z = 0 plane, two triangles per square.
    pts = zeros(n * n, 3);
    for j = 0:(n - 1)
        for i = 0:(n - 1)
            pts(j * n + i + 1, :) = [i, j, 0];
        end
    end
    conn = zeros(0, 3);
    for j = 0:(n - 2)
        for i = 0:(n - 2)
            a = j * n + i + 1;
            conn = [conn; a, a + 1, a + n + 1; a, a + n + 1, a + n]; %#ok<AGROW>
        end
    end
    m = mio.mesh();
    mio.set_points(m, pts);
    mio.add_cell_block(m, 'triangle', conn);
end

function m = sphere_cloud(n)
% A Fibonacci point cloud of N points on the unit sphere (no cells).
    k = (0:(n - 1))';
    z = 1 - (2 * k + 1) / n;
    r = sqrt(1 - z .^ 2);
    phi = k * pi * (3 - sqrt(5));
    m = mio.mesh();
    mio.set_points(m, [r .* cos(phi), r .* sin(phi), z]);
end

% --- tests ------------------------------------------------------------------

function t_geometry()
    m = fixture();
    surf = mio.extract_surface(m);
    assert(mio.num_points(surf) > 0);
    expect_equal(mio.cell_block_type(surf, 1), 'triangle');
    surf2 = mio.extract_surface(m, 'record_parent_ids', true);
    assert(any(strcmp(mio.cell_data_names(surf2), 'surface:parent_cell')));
    skin = mio.extract_skin(m, 'linearize', true);
    assert(mio.num_cells(skin) > 0);
    q = mio.attach_quality(m);
    assert(any(strncmp(mio.cell_data_names(q), 'quality:', 8)));
    c = mio.quality_counts(m);
    expect_equal(c.num_cells, 2);
    expect_equal(c.num_inverted, 0);
end

function t_transform()
    m = fixture();
    M = [1 0 0 2; 0 1 0 3; 0 0 1 4; 0 0 0 1];
    t = mio.transform(m, M);
    p = mio.points(t);
    expect_equal(p(1, :), [2 3 4]); % the origin moved by the translation
    fp = fixture_points();
    expect_equal(p(2, :), fp(2, :) + [2 3 4]);
    expect_equal(mio.connectivity(t, 1), fixture_conn()); % topology untouched
    % A non-symmetric linear block proves the row-/column-major transpose.
    R = [0 -1 0 0; 1 0 0 0; 0 0 1 0; 0 0 0 1];
    r = mio.transform(m, R, true);
    pr = mio.points(r);
    expect_equal(pr(2, :), [-0.2, 1.1, 0.3], 1e-12);
    expect_error(@() mio.transform(m, eye(3)), '4x4');
end

function t_clean_smooth_crop()
    m = fixture();
    c1 = mio.clean(m);
    expect_equal(mio.num_points(c1.mesh), 5); % nothing to remove here
    expect_equal(c1.points_welded, 0);

    s = mio.smooth(m, 'iterations', 2);
    expect_equal(mio.num_points(s.mesh), 5); % a pure coordinate move
    expect_equal(mio.num_cells(s.mesh), 2);
    assert(s.max_displacement >= 0);

    % method='odt': tet-only, C++-core only. The fixture's block is a tetra.
    so = mio.smooth(m, 'method', 'odt', 'iterations', 1, 'fix_boundary', false, ...
                    'preserve_features', false);
    expect_equal(mio.num_points(so.mesh), 5);
    expect_equal(mio.num_cells(so.mesh), 2);
    assert(so.num_nodes_moved >= 0);

    % optimize_volume: ODT remeshing (relocate + flip connectivity), tet-only.
    ov = mio.optimize_volume(m, 'max_iterations', 5);
    expect_equal(mio.num_points(ov.mesh), 5); % the point set is invariant
    assert(ov.num_flips >= 0);
    assert(ov.min_dihedral_after >= ov.min_dihedral_before);
    ov = mio.optimize_volume(m, 'max_iterations', 5, 'relocate', false, 'objective', 'energy');
    assert(ov.min_quality_after >= ov.min_quality_before - 1e-12);
    expect_equal(ov.num_edge_removals, 0);
    expect_error(@() mio.optimize_volume(m, 'objective', 'best'), 'objective');

    b = mio.crop_bbox(m, [-1 -1 -1], [10 10 10]);
    expect_equal(mio.num_cells(b), 2); % the box holds everything

    % frozen: reaches the core independently of fix_boundary.
    sf = mio.smooth(m, 'iterations', 3, 'fix_boundary', false, 'frozen', [1 2]);
    expect_equal(mio.num_points(sf.mesh), 5);
    % An out-of-range frozen id fails by name (0-based internally, so 1-based
    % id 999 is reported as 998).
    expect_error(@() mio.smooth(m, 'frozen', 999), 'frozen node id 998');

    b2 = mio.crop_bbox(m, [100 100 100], [101 101 101]);
    expect_equal(mio.num_cells(b2), 0);
    p = mio.crop_plane(m, [0 0 0], [1 0 0], 'mode', 'any');
    assert(mio.num_cells(p) >= 0);
    expect_error(@() mio.crop_bbox(m, [0 0 0], [1 1 1], 'mode', 'some'), 'mode');
end

function t_slice_iso()
    m = fixture();
    sl = mio.slice(m, [0.5 0.5 0.5], [0 0 1]);
    assert(mio.num_points(sl) >= 0); % may be empty; never an error
    iso = mio.isosurface(m, 'temperature', 3);
    assert(mio.num_points(iso) >= 0);
    if mio.num_cells(iso) > 0
        assert(any(strcmp(mio.cell_data_names(iso), 'iso:value')));
        assert(any(strcmp(mio.cell_data_names(iso), 'iso:index')));
    end
    % A cell_data array has no level set: that must fail by name.
    expect_error(@() mio.isosurface(m, 'material', 8));
    expect_error(@() mio.isosurface(m, 'temperature', []), 'isovalue');
end

function t_gradient()
    m = fixture();
    g = mio.gradient(m, 'temperature');
    expect_equal(g.num_skipped, 0);
    expect_equal(g.num_fallback, 0);
    assert(any(strcmp(mio.cell_data_names(g.mesh), 'temperature:gradient')));
    l = mio.gradient(m, 'temperature', 'method', 'least-squares');
    assert(l.num_fallback >= 0);
    p = mio.gradient(m, 'temperature', 'location', 'point', 'output', 'dT');
    assert(any(strcmp(mio.point_data_names(p.mesh), 'dT')));
    expect_equal(size(mio.point_data(p.mesh, 'dT')), [5 3]);
    % A cell_data array is piecewise constant and has no derivative; a scalar
    % has no divergence. Both must fail by name.
    expect_error(@() mio.gradient(m, 'material'));
    expect_error(@() mio.gradient(m, 'temperature', 'op', 'divergence'));
end

function t_hessian()
    m = fixture();
    h = mio.hessian(m, 'temperature');
    expect_equal(h.num_skipped, 0);
    assert(any(strcmp(mio.cell_data_names(h.mesh), 'temperature:hessian')));
    l = mio.hessian(m, 'temperature', 'method', 'least-squares');
    assert(l.num_fallback >= 0);
    p = mio.hessian(m, 'temperature', 'location', 'point', 'output', 'H2');
    assert(any(strcmp(mio.point_data_names(p.mesh), 'H2')));
    % A cell_data array has no derivative, and a vector field is scalar-only.
    expect_error(@() mio.hessian(m, 'material'));
    expect_error(@() mio.hessian(m, 'displacement'));
end

function t_estimate_error()
    m = fixture();
    e = mio.estimate_error(m, 'temperature');
    expect_equal(e.num_skipped, 0);
    assert(e.global_error >= 0);
    expect_equal(e.num_marked, 0);
    assert(any(strcmp(mio.cell_data_names(e.mesh), 'error:zz')));
    assert(~any(strcmp(mio.cell_data_names(e.mesh), 'error:marked')));
    f = mio.estimate_error(m, 'temperature', 'marking', 'absolute', 'marking_value', 1e-9, ...
                           'output', 'ind', 'marked', 'flag');
    assert(any(strcmp(mio.cell_data_names(f.mesh), 'ind')));
    assert(any(strcmp(mio.cell_data_names(f.mesh), 'flag')));
    assert(f.num_marked >= 0);
    % A cell_data array has no derivative to recover, and an out-of-range
    % marking_value for 'fraction' is rejected.
    expect_error(@() mio.estimate_error(m, 'material'));
    expect_error(@() mio.estimate_error(m, 'temperature', 'marking', 'fraction', ...
                                        'marking_value', 1.5));
end

function t_crop_predicate()
    % The R test crops a distance field on a grid (other files' functions);
    % the fixture's cell data does the same job here.
    m = fixture();
    kept = mio.crop_predicate(m, 'material', 'compare', '<', 'value', 8);
    info = mio.cell_block_info(kept, 1);
    expect_equal(info.num_cells, 1);
    expect_equal(mio.cell_data(kept, 'material'), 7);
    ne = mio.crop_predicate(m, 'material', [char(33) '='], 7);
    expect_equal(mio.cell_data(ne, 'material'), 9);
    all_kept = mio.crop_predicate(m, 'material', 'compare', '>=', 'value', 0);
    expect_equal(mio.num_cells(all_kept), 2);
    % There is deliberately no `mode`: a cell_data predicate is already one
    % value per cell and has nothing for an all/any rule to reduce.
    expect_error(@() mio.crop_predicate(m, 'material', 'compare', '~'), 'compare');
    expect_error(@() mio.crop_predicate(m, 'nope'));
end

function t_remesh()
    cube = cube_surface();
    r = mio.remesh(cube, 10);
    expect_equal(mio.cell_block_type(r.mesh, 1), 'triangle');
    expect_equal(mio.num_points(r.mesh), r.num_clusters);
    assert(r.num_isolated_clusters >= 0);
    assert(r.num_non_manifold_vertices >= 0);
    % gradation/preserve_boundary are accepted and reach the C API; a closed
    % mesh's boundary handling is a no-op either way.
    q = mio.remesh(cube, 10, 'metric', 'quadric', 'gradation', 1.5, 'preserve_boundary', false);
    assert(mio.num_cells(q.mesh) > 0);
    % The anisotropic metric + max_anisotropy go through mio_remesh_ex.
    a = mio.remesh(cube, 10, 'metric', 'anisotropic', 'max_anisotropy', 3.0);
    assert(mio.num_cells(a.mesh) > 0);
    assert(a.num_non_manifold_vertices >= 0);
    % max_anisotropy away from the default under a non-anisotropic metric is
    % guarded, not silently ignored.
    expect_error(@() mio.remesh(cube, 10, 'metric', 'isotropic', 'max_anisotropy', 2.0));
    expect_error(@() mio.remesh(cube, 3));
end

function t_normals()
    % The unit cube surface (12 consistently wound triangles): one smooth
    % normal per point, or 24 points once split at a 30 degree crease.
    m = cube_surface();
    s = mio.compute_normals(m);
    expect_equal(mio.num_points(s.mesh), 8);
    expect_equal(s.num_added_points, 0);
    assert(s.quality.watertight);
    expect_equal(size(mio.point_data(s.mesh, 'normals')), [8 3]);
    sp = mio.compute_normals(m, 'split_angle', 30, 'cell_normals', true, ...
                             'record_parent_ids', true);
    expect_equal(mio.num_points(sp.mesh), 24);
    expect_equal(sp.num_added_points, 16);
    expect_equal(sp.num_split_points, 8);
    expect_error(@() mio.compute_normals(m, 'weight', 'bogus'), 'unknown weight');
    expect_error(@() mio.compute_normals(m, 'split_angle', 270));
end

function t_curvature()
    m = cube_surface();
    c = mio.compute_curvature(m, 'record_area', true, 'record_principal', true);
    % On a closed surface the total angle defect is 2*pi*chi = 4*pi.
    expect_equal(c.total_angle_defect, 4 * pi, 1e-9);
    assert(c.quality.watertight);
    names = mio.point_data_names(c.mesh);
    assert(any(strcmp(names, 'curvature:mean')));
    assert(any(strcmp(names, 'curvature:gaussian')));
    assert(any(strcmp(names, 'curvature:area')));
    expect_equal(size(mio.point_data(c.mesh, 'curvature:principal')), [8 2]);
    b = mio.compute_curvature(m, 'dual_area', 'barycentric');
    expect_equal(b.total_angle_defect, 4 * pi, 1e-9);
    expect_error(@() mio.compute_curvature(m, 'dual_area', 'nope'), 'unknown dual area');
end

function t_repair()
    % The cube surface with two facets flipped: exactly those two are rewound
    % and the output is watertight with no point added.
    cube = cube_surface();
    pts = mio.points(cube);
    conn = mio.connectivity(cube, 1);
    broken = conn;
    broken(1, [2 3]) = broken(1, [3 2]);
    broken(9, [2 3]) = broken(9, [3 2]);
    m = mio.mesh();
    mio.set_points(m, pts);
    mio.add_cell_block(m, 'triangle', broken);
    r = mio.repair(m);
    assert(r.quality_before.inconsistent_pairs > 0);
    assert(r.quality_after.watertight);
    expect_equal(r.num_flipped, 2);
    expect_equal(r.num_components, 1);
    expect_equal(r.num_holes_detected, 0);
    expect_equal(mio.num_points(r.mesh), 8);

    % Minus one facet: the hole is filled by a centroid fan in a trailing block.
    h = mio.mesh();
    mio.set_points(h, pts);
    mio.add_cell_block(h, 'triangle', conn(2:end, :));
    f = mio.repair(h);
    expect_equal(f.num_holes_filled, 1);
    expect_equal(f.num_faces_added, 3);
    expect_equal(f.num_points_added, 1);
    assert(f.quality_after.watertight);
    expect_equal(mio.num_cell_blocks(f.mesh), 2);
    g = mio.repair(h, 'fill_holes', false);
    expect_equal(g.num_holes_detected, 0);
    expect_equal(g.quality_after.boundary_edges, 3);
end

function t_shrinkwrap()
    target = mio.mesh();
    mio.set_points(target, [0 0 0; 4 0 0; 4 4 0; 0 4 0]);
    mio.add_cell_block(target, 'triangle', [1 2 3; 1 3 4]);
    cloud = mio.mesh();
    mio.set_points(cloud, [1 1 0.5; 2 2 -0.7; 3 1 1.5]);
    mio.add_cell_block(cloud, 'vertex', [1; 2; 3]);
    mio.add_point_data(cloud, 'w', [1; 0.5; 1]);
    w = mio.shrinkwrap(cloud, target, 'offset', 0.25, 'weights', 'w', 'record_distance', true);
    expect_equal(w.num_projected, 3);
    expect_equal(w.num_skipped, 0);
    expect_equal(w.quality.boundary_edges, 4);
    p = mio.points(w.mesh);
    expect_equal(p(1, 3), 0.25, 1e-14);
    expect_equal(p(2, 3), -0.7 + 0.5 * (0.25 + 0.7), 1e-14);
    assert(any(strcmp(mio.point_data_names(w.mesh), 'shrinkwrap:distance')));
    expect_error(@() mio.shrinkwrap(cloud, target, 'normal_weight', 'nope'));
end

function t_sobolev()
    n = 4;
    m = tri_grid(n);
    pts = mio.points(m);
    disp = zeros(n * n, 3);
    for j = 0:(n - 1)
        for i = 0:(n - 1)
            if mod(i + j, 2) == 1
                disp(j * n + i + 1, 3) = 0.1;
            else
                disp(j * n + i + 1, 3) = -0.1;
            end
        end
    end
    mio.add_point_data(m, 'd', disp);
    pin = zeros(n * n, 1);
    pin(1) = 1;
    mio.add_point_data(m, 'pin', pin);
    s = mio.sobolev_deform(m, 'd', 2.0, 'fixed_points_array', 'pin', 'record_filtered', true);
    assert(s.converged);
    assert(s.num_iterations > 0);
    expect_equal(s.num_fixed, 1);
    expect_equal(s.num_isolated, 0);
    assert(s.max_displacement > 0 && s.max_displacement < 0.1);
    ps = mio.points(s.mesh);
    expect_equal(ps(1, 3), 0);
    assert(any(strcmp(mio.point_data_names(s.mesh), 'sobolev:displacement')));
    % A constant field with nothing pinned is preserved exactly, zero
    % iterations.
    mio.add_point_data(m, 'c', 0.4 * ones(n * n, 3));
    k = mio.sobolev_deform(m, 'c', 2.0);
    expect_equal(k.num_iterations, 0);
    pk = mio.points(k.mesh);
    expect_equal(pk(:, 1), pts(:, 1) + 0.4);
    expect_error(@() mio.sobolev_deform(m, 'missing', 1.0));
end

function t_feature_segment()
    m = cube_surface();
    fe = mio.feature_edges(m);
    expect_equal(fe.num_feature, 12);
    expect_equal(fe.num_boundary, 0);
    expect_equal(mio.cell_block_type(fe.mesh, 1), 'line');
    expect_error(@() mio.feature_edges(m, 'feature_angle', 200));

    sg = mio.segment(m);
    expect_equal(sg.num_patches, 6);
    fp = mio.fit_primitives(sg.mesh);
    expect_equal(fp.num_plane, 6);
    mz = mio.fit_primitives(sg.mesh, 'criterion', 'minimum_zone');
    expect_equal(mz.num_plane, 6);
    expect_error(@() mio.fit_primitives(sg.mesh, 'criterion', 'chebyshev'), 'criterion');
    step_path = [tempname() '.step'];
    el = mio.write_step_elements(step_path, fp.mesh);
    expect_equal(el.num_planes, 6);
    el = mio.write_step_elements(step_path, fp.mesh, 'patches', 0);
    expect_equal(el.num_planes, 1);
    expect_error(@() mio.write_step_elements(step_path, fp.mesh, 'patches', 99), 'patch');
    cl = mio.clip(m, [0 0 0.5], [0 0 1]);
    expect_equal(cl.num_loops, 1);
    fl = mio.flip_normals(m);
    expect_equal(fl.num_flipped, 12);
    ri = mio.remove_islands(m);
    expect_equal(ri.num_components, 1);
    expect_error(@() mio.remove_islands(m, 'connectivity', 'face'), 'connectivity');
    ms = mio.measure(m, ['{"elements": [' ...
        '{"name": "A", "fit": "plane", "seeds": [[0.5, 0.5, 0.0]]},' ...
        '{"name": "B", "fit": "plane", "seeds": [[1.0, 0.5, 0.5]]}],' ...
        '"tolerances": [{"type": "perpendicularity", "feature": "B", "datum": "A", "limit": 1e-9}]}']);
    expect_equal(ms.num_elements, 2);
    expect_equal(ms.num_fail, 0);
    al = mio.align_datums(ms.mesh, 'A', 'B');
    expect_equal(size(al.matrix), [4 4]);
    el = mio.measure_element_primitive(ms.mesh, 'A');
    expect_equal(el.type, 1);
    delete(step_path);
    fp2 = mio.fit_primitives(sg.mesh, 'types', 'plane');
    expect_equal(fp2.num_plane, 6);
    expect_error(@() mio.segment(m, 'feature_angle', 200));
    expect_error(@() mio.segment(m, 'method', 'blob'), 'method');
    expect_error(@() mio.fit_primitives(m));
    expect_error(@() mio.fit_primitives(sg.mesh, 'types', 'blob'), 'blob');
    fp3 = mio.fit_primitives(sg.mesh, 'types', 'plane', 'seeds', [0.5 0.5 1]);
    expect_equal(fp3.num_plane, 1);
    assert(fp3.num_feature_cells > 0);
    expect_error(@() mio.fit_primitives(sg.mesh, 'used_points', '5sigma'), 'used_points');
end

function t_refine_patches()
    m = cube_surface();
    sg = mio.segment(m);
    % refine_patches (v17.11.0): six planes are already one patch each.
    rp = mio.refine_patches(sg.mesh);
    expect_equal(rp.num_patches, 6);
    expect_equal(rp.num_primitive, 6);
    expect_equal(rp.num_splits + rp.num_merges, 0);
    expect_equal(rp.explained_area_fraction, 1, 1e-12);
    assert(any(strcmp(rp.stop, {'unchanged', 'small_gain', 'iteration_cap', 'energy_rose'})));
    rp0 = mio.refine_patches(sg.mesh, 'max_iterations', 0, 'ransac_seed', 42, 'types', 'plane');
    expect_equal(rp0.stop, 'iteration_cap');
    expect_equal(rp0.num_iterations, 0);
    expect_error(@() mio.refine_patches(m), 'segment');
    expect_error(@() mio.refine_patches(sg.mesh, 'tolerance', 0));
    expect_error(@() mio.refine_patches(sg.mesh, 'types', 'blob'));
    expect_error(@() mio.refine_patches(sg.mesh, 'ransac_seed', -1), 'ransac_seed');
    expect_error(@() mio.refine_patches(sg.mesh, 'ransac_seed', 1.5), 'ransac_seed');
end

function t_fit_cad()
    m = cube_surface();
    assert(islogical(mio.has_occt()));
    % fit_cad refines the segmentation on request, on every backend.
    rc = mio.fit_cad(m, 'backend', 'native', 'refine_segmentation', true, ...
                     'max_refine_iterations', 5);
    assert(rc.refined);
    expect_equal(rc.num_plane, 6);
    expect_equal(rc.explained_area_fraction, 1, 1e-12);
    pc = mio.fit_cad(m, 'backend', 'native');
    assert(~pc.refined);
    assert(~pc.refinement_kept);
    expect_equal(pc.num_refine_splits, 0);
    expect_equal(pc.backend, 'native');
    expect_equal(pc.freeform_mode, 'legacy');
    % The Phase 3 report (v17.11.0): closed on every backend, NaN unmeasured.
    assert(pc.closed);
    assert(isnan(pc.hausdorff_max));
    expect_equal(pc.num_snapped_junctions, 0);
    expect_error(@() mio.fit_cad(m, 'backend', 'native', 'max_freeform_poles', 3), ...
                 'max_freeform_poles');
    expect_error(@() mio.fit_cad(m, 'backend', 'native', 'sew_tolerance', -1), 'sew_tolerance');
    gc = mio.fit_cad(m, 'backend', 'native', 'freeform_general', true, 'max_freeform_splits', 0);
    assert(gc.closed);
    expect_equal(gc.num_freeform_general, 0);
    % The constrained primitives: off by default, the cube keeps its faces.
    expect_equal(pc.num_constrained_primitives, 0);
    expect_equal(pc.num_constraint_groups, 0);
    cc = mio.fit_cad(m, 'backend', 'native', 'constrain_primitives', true);
    expect_equal(cc.num_faces, 6);
    expect_equal(cc.num_demoted, 0);
    assert(cc.num_constrained_primitives >= 0 && cc.num_constrained_primitives <= 6);
    % A STEP file is written on request.
    step = [tempname() '.step'];
    cleanup = onCleanup(@() delete_if(step));
    st = mio.fit_cad(m, 'step_path', step, 'backend', 'native');
    expect_equal(st.num_plane, 6);
    assert(exist(step, 'file') == 2);
    if mio.has_occt()
        sc = mio.fit_cad(m, 'snap_junctions', true, 'face_deviation', true);
        assert(sc.closed);
        assert(sc.num_snapped_junctions > 0);
        assert(isfinite(sc.hausdorff_max));
        cad = mio.fit_cad(m);
        assert(cad.valid);
        expect_equal(cad.num_plane, 6);
        expect_equal(cad.num_split_patches, 0);
        expect_equal(cad.num_fillet, 0);
        expect_equal(cad.num_helical, 0);
        expect_equal(cad.num_variable_fillet, 0);
        expect_equal(cad.num_blend, 0);
        expect_equal(cad.num_rounded_thread, 0);
        expect_equal(cad.num_vertex_blend, 0);
        expect_equal(cad.num_chamfer, 0);
        expect_equal(cad.num_conic_blend, 0);
        expect_equal(cad.num_tapered_thread, 0);
        expect_equal(cad.num_variable_pitch_thread, 0);
        bare = mio.fit_cad(m, 'freeform_split', false, 'fillets', false, 'threads', false);
        expect_equal(bare.num_plane, 6);
    else
        fallback = mio.fit_cad(m);
        expect_equal(fallback.backend, 'native');
        expect_error(@() mio.fit_cad(m, 'backend', 'occt'), 'MESHIOPLUSPLUS_WITH_OCCT');
    end
    expect_error(@() mio.fit_cad(m, 'types', 'blob'));
    expect_error(@() mio.fit_cad(m, 'freeform_mode', 'nope'), 'freeform_mode');
    expect_error(@() mio.fit_cad(m, 'backend', 'nope'), 'backend');
end

function t_hausdorff_quality()
    m = cube_surface();
    h = mio.hausdorff_distance(m, m, 'face_samples', 2);
    expect_equal(h.distance, 0);
    expect_equal(h.num_samples_a, 8 + 12 * 4);
    expect_equal(size(h.worst_point_a), [1 3]);
    moved = mio.transform(m, [1 0 0 0; 0 1 0 0; 0 0 1 0.5; 0 0 0 1]);
    d = mio.hausdorff_distance(m, moved);
    expect_equal(d.distance, 0.5, 1e-12);

    q = mio.check_quality(m, 'min_angle >= 30');
    assert(q.passed);
    expect_equal(q.num_checks, 3);
    assert(ischar(q.summary) && ~isempty(q.summary));
    q2 = mio.check_quality(m, 'min_angle >= 50');
    assert(~q2.passed);
    expect_error(@() mio.check_quality(m, 'bogus >= 1'));
end

function t_geodesic_etc()
    m = cube_surface();
    % Geodesic distance from one cube corner: the far one is about 2.24 away.
    gd = mio.geodesic_distance(m, 1);
    assert(gd.converged);
    expect_equal(gd.num_unreached, 0);
    assert(gd.max_distance > 1.8 && gd.max_distance < 2.6);
    expect_equal(numel(mio.point_data(gd.mesh, 'geodesic:distance')), 8);
    expect_error(@() mio.geodesic_distance(m), 'source');
    mio.add_region(m, 'corner', 'point', 1);
    gr = mio.geodesic_distance(m, 'source_region', 'corner');
    expect_equal(gr.num_sources, 1);

    % A UV atlas of the closed cube: several charts, nothing flipped.
    pa = mio.parameterize(m);
    assert(pa.num_charts >= 2);
    expect_equal(pa.num_flipped, 0);
    assert(pa.num_seam_points > 0);
    expect_error(@() mio.parameterize(m, 'max_stretch', 0.5), 'max_stretch');

    expect_error(@() mio.morph(m, 'missing'));
    expect_error(@() mio.untangle(m));
    expect_error(@() mio.prepare_metric(m, 'missing'));

    % Smoke: a zero displacement morphs a tet mesh to itself.
    t = fixture();
    mio.add_point_data(t, 'zero', zeros(5, 3));
    mo = mio.morph(t, 'zero');
    expect_equal(mio.points(mo.mesh), fixture_points(), 1e-12);
    expect_equal(mo.num_inverted, 0);
    % Smoke: an untangled tet mesh stays untangled.
    un = mio.untangle(t);
    expect_equal(un.num_inverted_after, 0);
    % Smoke: an isotropic metric field (6 components) prepares.
    mio.add_point_data(t, 'met', repmat([1 0 0 1 0 1], 5, 1));
    pm = mio.prepare_metric(t, {'met'}, 'output', 'M');
    assert(any(strcmp(mio.point_data_names(pm.mesh), 'M')));
    assert(isfield(pm, 'complexity_after') && pm.num_graded >= 0);
end

function t_spectral_registration()
    m = cube_surface();
    ee = mio.laplacian_eigenmodes(m, 'num_modes', 4);
    assert(ee.converged);
    expect_equal(size(mio.point_data(ee.mesh, 'laplacian:modes')), [8 4]);
    expect_error(@() mio.laplacian_eigenmodes(m, 'num_modes', 9));

    rr = mio.register_rigid(m, m);
    assert(rr.converged);
    expect_equal(rr.matrix, eye(4), 1e-12);
    expect_error(@() mio.register_rigid(m, m, 'method', 'affine'), 'method');
    % A translated copy registers back; the matrix is a mathematical 4x4, so
    % the translation sits in its last column.
    moved = mio.transform(m, [1 0 0 0.05; 0 1 0 0; 0 0 1 0; 0 0 0 1]);
    rt = mio.register_rigid(moved, m);
    expect_equal(rt.matrix(1, 4), -0.05, 1e-6);
    rs = mio.register_rigid(moved, m, 'initial', eye(4));
    expect_equal(rs.matrix(1, 4), -0.05, 1e-6);
    expect_error(@() mio.register_rigid(m, m, 'initial', eye(3)), '4 x 4');
    % The automatic start, landmarks, the symmetry plane and auto-align
    % (roadmap 7.3).
    ra = mio.register_rigid(moved, m, 'initial', 'auto', 'symmetric_check', true);
    assert(ra.rms < 1e-6);
    expect_equal(numel(ra.candidate_rms), 24);
    expect_error(@() mio.register_rigid(m, m, 'initial', 'guess'), 'initial');
    lm = mio.register_landmarks([0 0 0; 1 0 0; 0 2 0], [1 2 3; 2 2 3; 1 4 3]);
    expect_equal(lm.matrix(1:3, 4), [1; 2; 3], 1e-12);
    expect_error(@() mio.register_landmarks([0 0 0; 1 0 0], [0 0 0; 1 0 0]));
    sp = mio.find_symmetry_plane(m);
    assert(sp.share > 0.9);
    al = mio.auto_align(m);
    assert(strcmp(al.mode, 'face'));
    dv = mio.deviation(moved, 'reference', m, 'register', true);
    assert(dv.registered);

    nn = mio.estimate_normals(m, 'k', 3);
    expect_equal(nn.num_degenerate, 0);
    expect_equal(size(mio.point_data(nn.mesh, 'normals')), [8 3]);
    expect_error(@() mio.estimate_normals(m, 'orientation', 'inward'), 'orientation');
    expect_error(@() mio.estimate_normals(m, 'orientation', 'viewpoint', 'viewpoint', [1 2]), ...
                 'three numbers');
end

function t_reconstruct()
    cloud = sphere_cloud(800);
    rs = mio.reconstruct_surface(cloud);
    assert(rs.estimated_normals);
    expect_equal(rs.num_components, 1);
    assert(rs.num_triangles > 0);
    expect_error(@() mio.reconstruct_surface(cloud, 'k', 0), 'k must be');
    expect_error(@() mio.reconstruct_surface(cloud, 'estimate_normals', false), 'normals');
end

function t_hull_coords()
    m = cube_surface();
    hh = mio.convex_hull(m);
    expect_equal(hh.num_vertices, 8);
    expect_equal(hh.volume, 1, 1e-12);
    expect_equal(hh.dimension, 3);
    dd = mio.delaunay_2d(m);
    expect_equal(dd.num_triangles, 2);
    expect_error(@() mio.convex_hull(m, 5), 'dimension');
    cc = mio.convert_coordinates(m, 'spherical', 'origin', [0.5 0.5 0.5]);
    expect_equal(cc.num_singular, 0);
    expect_equal(size(mio.point_data(cc.mesh, 'coords:spherical')), [8 3]);
    expect_error(@() mio.convert_coordinates(m, 'polar'), 'cylindrical');
    expect_error(@() mio.convert_coordinates(m, 'axis', [0 1]), 'three numbers');
    % A vector array is rewritten in the local basis.
    mio.add_point_data(m, 'v', repmat([1 0 0], 8, 1));
    cv = mio.convert_coordinates(m, 'arrays', {'v'});
    expect_equal(cv.num_arrays, 1);
end

function t_interface()
    m = mio.mesh();
    mio.set_points(m, [0 0 0; 1 0 0; 0 1 0; 0 0 1; 0 0 -1]);
    mio.add_cell_block(m, 'tetra', [1 2 3 4; 1 3 2 5]);
    mio.add_region(m, 'upper', 'cell', 1, 'dim', 3);
    mio.add_region(m, 'lower', 'cell', 2, 'dim', 3);
    mio.add_region(m, 'slave', 'point', 1);
    mio.add_region(m, 'cut', 'side', [1 3]);

    found = mio.find_interface(m, 'upper', 'lower');
    expect_equal(found.num_pairs, 1);
    expect_equal(found.side_a, [1 3]); % 1-based cell, 0-based facet
    expect_equal(size(found.side_b), [1 2]);
    expect_equal(mio.cell_block_type(found.mesh, 1), 'triangle');
    expect_error(@() mio.find_interface(m, 'upper', 'lower', 'mode', 'nope'), 'mode');
    expect_error(@() mio.find_interface(m, 'upper', 'lower', 'master', 'c'), 'master');

    contacts = mio.contact_pairs(m, 'slave', 'lower', 'tolerance', 10);
    expect_equal(contacts.slave_point, 1);
    expect_equal(contacts.master_cell, 2);
    expect_equal(size(contacts.local_coordinates), [1 3]);
    expect_equal(size(contacts.normal), [1 3]);
    expect_equal(numel(contacts.unmatched), 0);

    sp = mio.split_interface(m, 'cut', 'add_cohesive', true);
    expect_equal(sp.num_duplicated_points, 3);
    expect_equal(sp.num_cohesive_cells, 1);
end

function t_region_adjacency()
    m = mio.mesh();
    mio.set_points(m, [0 0 0; 1 0 0; 0 1 0; 0 0 1; 0 0 -1]);
    mio.add_cell_block(m, 'tetra', [1 2 3 4; 1 3 2 5]);
    mio.add_region(m, 'upper', 'cell', 1, 'dim', 3);
    mio.add_region(m, 'lower', 'cell', 2, 'dim', 3);
    a = mio.region_adjacency(m);
    expect_equal(mio.num_cells(a), 1);
    b = mio.region_adjacency(m, {'upper', 'lower'});
    expect_equal(mio.num_cells(b), 1);
end

function delete_if(path)
    if exist(path, 'file')
        delete(path);
    end
end
