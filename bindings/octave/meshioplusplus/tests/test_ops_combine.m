function tests = test_ops_combine()
%TEST_OPS_COMBINE  Combining, comparing, statistics and the opaque-result
%   operations (ops_combine.c): merge, interpolate, diff, stats, split,
%   partition, reorder, refine, decimate, sweeps, MMG and the rest.
    tests = {
        'stats and bandwidth summarise the fixture', @t_stats
        'tessellate splits a tetra10 conformingly', @t_tessellate
        'refine reports 1-based maps with a 0 sentinel', @t_refine_maps
        'subdivide polyhedrally refines with no point map', @t_subdivide
        'tet_to_poly builds the dual and poly_to_tet decomposes it', @t_tet_poly
        'agglomerate coarsens with a flat cell map', @t_agglomerate
        'selective refine closes up conformingly', @t_refine_selective
        'refine hierarchy is opt-in and maintained', @t_refine_hierarchy
        'undo_green restores the coarse parent', @t_undo_green
        'split and partition pieces own their handles', @t_split_partition
        'reorder returns 1-based permutations', @t_reorder
        'decimate honours target_faces and frozen', @t_decimate
        'decimate_volume collapses tets and reports its maps', @t_decimate_volume
        'error_size turns the indicator into a size field', @t_error_size
        'check_tetrahedra counts the tetrahedra', @t_check_tetrahedra
        'check_conformity reports a conforming mesh', @t_conformity
        'merge and interpolate work', @t_merge_interpolate
        'diff reports structured results', @t_diff
        'blend_steps blends the data of two steps', @t_blend
        'remesh_quads lists and names its backends', @t_remesh_quads
        'cut_level_set and cut_quadrature cut the cube skin', @t_cut
        'boundary_layer grows wedges from a surface', @t_boundary_layer
        'extrude, revolve and mirror sweep and reflect', @t_sweeps
        'MMG operations run or name the build flag', @t_mmg
        'metric_from_hessian writes a metric', @t_metric
        'tetrahedralize fills a closed surface', @t_tetrahedralize
        'find_intersections and thickness', @t_intersections
        'deviation measures an offset sphere', @t_deviation
        'capability queries return their types', @t_capabilities
    };
end

% --- local fixtures --------------------------------------------------------

function tf = has_cmd(name)
    tf = any(strcmp(meshioplusplus_mex('__commands__'), name));
end

function m = grid2d(nx, ny, cell_type)
%   The unit square at z = 0 as nx x ny quads (or 2 triangles each), points
%   numbered x fastest -- the shape of mio.rectangle(dims = [nx ny]).
    [X, Y] = meshgrid((0:nx) / nx, (0:ny) / ny);
    X = X.';
    Y = Y.';
    p = [X(:), Y(:), zeros(numel(X), 1)];
    id = @(i, j) j * (nx + 1) + i + 1;
    quads = zeros(nx * ny, 4);
    k = 0;
    for j = 0:ny - 1
        for i = 0:nx - 1
            k = k + 1;
            quads(k, :) = [id(i, j), id(i + 1, j), id(i + 1, j + 1), id(i, j + 1)];
        end
    end
    m = mio.mesh();
    mio.set_points(m, p);
    if strcmp(cell_type, 'triangle')
        mio.add_cell_block(m, 'triangle', [quads(:, [1 2 3]); quads(:, [1 3 4])]);
    else
        mio.add_cell_block(m, 'quad', quads);
    end
end

function m = hex_box(n)
%   The unit cube as n^3 hexahedra.
    [X, Y, Z] = ndgrid((0:n) / n, (0:n) / n, (0:n) / n);
    p = [X(:), Y(:), Z(:)];
    id = @(i, j, k) k * (n + 1)^2 + j * (n + 1) + i + 1;
    hex = zeros(n^3, 8);
    c = 0;
    for k = 0:n - 1
        for j = 0:n - 1
            for i = 0:n - 1
                c = c + 1;
                hex(c, :) = [id(i, j, k), id(i + 1, j, k), id(i + 1, j + 1, k), id(i, j + 1, k), ...
                             id(i, j, k + 1), id(i + 1, j, k + 1), id(i + 1, j + 1, k + 1), ...
                             id(i, j + 1, k + 1)];
            end
        end
    end
    m = mio.mesh();
    mio.set_points(m, p);
    mio.add_cell_block(m, 'hexahedron', hex);
end

function m = fan()
%   A 4-triangle fan around a centre vertex (point 5), the fixture the
%   Fortran and C API tests use.
    m = mio.mesh();
    mio.set_points(m, [0 0 0; 1 0 0; 1 1 0; 0 1 0; 0.5 0.5 0]);
    mio.add_cell_block(m, 'triangle', [1 2 5; 2 3 5; 3 4 5; 4 1 5]);
end

% --- tests -----------------------------------------------------------------

function t_stats()
    m = fixture();
    st = mio.stats(m);
    expect_equal(st.num_points, 5);
    expect_equal(st.num_cells, 2);
    expect_equal(size(st.bbox_min), [1, 3]);
    expect_equal(st.bbox_min(1), 0.0);
    expect_equal(st.bbox_max(1), 1.4);
    expect_equal(st.bbox_max, max(fixture_points(), [], 1));
    assert(mio.compute_bandwidth(m) >= 0);
end

function t_tessellate()
    m = mio.mesh();
    mio.set_points(m, [0 0 0; 1 0 0; 0 1 0; 0 0 1; 0.5 0 0; 0.5 0.5 0; 0 0.5 0; ...
                       0 0 0.5; 0.5 0 0.5; 0 0.5 0.5]);
    mio.add_cell_block(m, 'tetra10', 1:10);
    t = mio.tessellate(m);
    expect_equal(t.num_curved_source_cells, 1);
    expect_equal(t.num_points, 35);
    expect_equal(mio.num_points(t.mesh), 35);
    expect_equal(t.num_cells, 64);
    expect_equal(t.num_boundary_facets, 64);
    expect_equal(t.num_interior_facets, 96);
    expect_equal(t.num_facets_with_bad_count, 0);
    src = double(mio.point_data(t.mesh, 'tessellate:source_point'));
    expect_equal(src(1:5), [0; 1; 2; 3; -1]);
    s = mio.tessellate(m, 'levels', 1, 'record_stencil', true);
    expect_equal(s.num_cells, 8);
    expect_equal(size(mio.point_data(s.mesh, 'tessellate:weights')), [10, 27]);
    p = mio.tessellate(m, 'curved', false);
    expect_equal(p.num_pass_through_source_cells, 1);
    expect_error(@() mio.tessellate(m, 'levels', -1));
end

function t_refine_maps()
    m = fixture();
    r = mio.refine(m, 'levels', 1);
    expect_equal(mio.num_cells(r.mesh), 16); % 2 tetra -> 8 children each
    % Refinement never prunes, so the point map is the identity -- 1-based.
    expect_equal(r.point_map, (1:5).');
    expect_equal(numel(r.cell_maps), 1);
    expect_equal(r.cell_maps{1}, [1; 9]); % each parent's first child
    % convert_cells CAN prune, and a pruned point is reported as 0.
    cc = mio.convert_cells(m, 'simplexify');
    expect_equal(mio.num_cells(cc.mesh), 2); % a tetra is already a simplex
    assert(all(cc.point_map >= 0));
    expect_equal(numel(cc.cell_maps), 1);
end

function t_subdivide()
    m = fixture();
    s = mio.subdivide(m, 'record_parent_ids', true);
    expect_equal(mio.num_cell_blocks(s.mesh), 1);
    expect_equal(mio.cell_block_type(s.mesh, 1), 'polyhedron');
    info = mio.cell_block_info(s.mesh, 1);
    assert(info.is_polyhedron);
    expect_equal(info.num_cells, 8); % 2 tetra, 4 faces each
    expect_equal(mio.num_points(s.mesh), mio.num_points(m) + 2); % one apex each
    assert(~isfield(s, 'point_map'));
    expect_equal(numel(s.cell_maps), 1);
    expect_equal(s.cell_maps{1}, [1; 5]);
    assert(any(strcmp('subdivide:parent_cell', mio.cell_data_names(s.mesh))));
end

function t_tet_poly()
    m = fixture(); % 2 tetra sharing one face, 5 vertices
    d = mio.tet_to_poly(m, 'record_vertex_ids', true);
    expect_equal(mio.num_cell_blocks(d.mesh), 1);
    expect_equal(mio.cell_block_type(d.mesh, 1), 'polyhedron');
    info = mio.cell_block_info(d.mesh, 1);
    expect_equal(info.num_cells, 5); % one per vertex
    expect_equal(sort(d.vertex_map), (1:5).');
    assert(d.num_new_points > 0);
    expect_equal(d.num_split_vertices, 0);
    b = mio.poly_to_tet(d.mesh);
    expect_equal(mio.cell_block_type(b.mesh, 1), 'tetra');
    expect_equal(b.num_exact_tets, 0);
    expect_equal(b.num_fanned_cells, 5);
    assert(iscell(b.cell_maps));
end

function t_agglomerate()
    m = fixture(); % 2 tetra sharing one face
    a = mio.agglomerate(m, 'target_group_size', 2);
    expect_equal(mio.num_cell_blocks(a.mesh), 1);
    expect_equal(mio.cell_block_type(a.mesh, 1), 'polyhedron');
    info = mio.cell_block_info(a.mesh, 1);
    assert(info.is_polyhedron);
    expect_equal(info.num_cells, 1); % both tetra merge into one cell
    expect_equal(mio.num_points(a.mesh), mio.num_points(m)); % never pruned
    % cell_map is a single FLAT array: both input cells land in one output cell.
    expect_equal(numel(a.cell_map), 2);
    expect_equal(a.cell_map(1), a.cell_map(2));
    b = mio.agglomerate(m, 'target_group_size', 2, 'merge_coplanar_faces', true, ...
                        'min_sphericity', 0.1);
    assert(b.num_faces_merged >= 0);
    assert(b.num_rejected >= 0);
end

function t_refine_selective()
    m = fixture();
    % `cells` is 1-based here and shifted to 0-based inside.
    sel = mio.refine(m, 'cells', 1, 'record_levels', true);
    assert(mio.num_cells(sel.mesh) > mio.num_cells(m));
    assert(mio.num_cells(sel.mesh) < 16); % not the uniform 8-per-cell
    assert(any(strcmp('refine:level', mio.cell_data_names(sel.mesh))));
    % Propagation reaches every cell of a connected mesh: uniform refinement.
    prop = mio.refine(m, 'cells', 1, 'closure', 'propagate');
    expect_equal(mio.num_cells(prop.mesh), 16);
    % At most one selector, and both enum names are validated.
    expect_error(@() mio.refine(m, 'cells', 1, 'region', 'anything'));
    expect_error(@() mio.refine(m, 'closure', 'blue'), 'closure');
    expect_error(@() mio.refine(m, 'where_array', 'q', 'where_op', '~='), 'comparison');
end

function t_refine_hierarchy()
    m = fixture();
    plain = mio.refine(m, 'cells', 1);
    assert(~any(strcmp('refine:cell_id', mio.cell_data_names(plain.mesh))));
    hier = mio.refine(m, 'cells', 1, 'record_hierarchy', true);
    assert(any(strcmp('refine:cell_id', mio.cell_data_names(hier.mesh))));
    assert(any(strcmp('refine:parent_id', mio.cell_data_names(hier.mesh))));
    ids = double(mio.cell_data(hier.mesh, 'refine:cell_id', 1));
    parents = double(mio.cell_data(hier.mesh, 'refine:parent_id', 1));
    % Stable identifiers on the raw C numbering, not shifted index maps.
    expect_equal(numel(unique(ids)), numel(ids));
    assert(all(ismember(parents, [0 1])));
    assert(any(strcmp('refine:entity', mio.point_data_names(hier.mesh))));
    % An existing hierarchy is maintained without the flag.
    again = mio.refine(hier.mesh, 'cells', 1);
    assert(any(strcmp('refine:cell_id', mio.cell_data_names(again.mesh))));
end

function t_undo_green()
    m = fixture(); % 2 tetrahedra sharing a whole face
    fine = mio.refine(m, 'cells', 1, 'record_hierarchy', true, 'record_levels', true);
    i_fine = mio.cell_block_info(fine.mesh, 1);
    i_coarse = mio.cell_block_info(m, 1);
    undone = mio.undo_green(m, fine.mesh);
    assert(undone.num_groups_undone > 0);
    assert(undone.num_cells_removed > 0);
    i_undone = mio.cell_block_info(undone.mesh, 1);
    assert(i_undone.num_cells < i_fine.num_cells);
    assert(i_undone.num_cells > i_coarse.num_cells);
    assert(~any(strcmp('refine:cell_id', mio.cell_data_names(undone.mesh))));
    assert(~any(strcmp('refine:entity', mio.point_data_names(undone.mesh))));
    % Fails by name rather than guessing: no hierarchy at all here.
    expect_error(@() mio.undo_green(m, m));
end

function t_split_partition()
    m = fixture();
    [pieces, keys] = mio.split(m, 'by', 'type');
    expect_equal(size(pieces), [1, 1]);
    expect_equal(keys, {'tetra'});
    % The piece still works after the C result is gone: it owns its handle.
    expect_equal(mio.num_cells(pieces{1}), 2);
    expect_equal(mio.points(pieces{1}), fixture_points());
    mio.release(pieces{1});

    parts = mio.partition(m, 2, 'method', 'sfc');
    expect_equal(size(parts), [2, 1]);
    expect_equal([parts.part_id], [0, 1]);
    total = 0;
    for k = 1:numel(parts)
        total = total + mio.num_cells(parts(k).mesh);
        expect_equal(numel(parts(k).cell_maps), 1); % block structure kept 1:1
    end
    expect_equal(total, 2); % partition of unity

    labels = mio.partition_labels(m, 2, 'method', 'sfc');
    expect_equal(size(labels), [2, 1]);
    % Part IDs, not indices: deliberately NOT shifted to 1-based.
    expect_equal(sort(labels), [0; 1]);
end

function t_reorder()
    m = fixture();
    r = mio.reorder(m, 'rcm');
    expect_equal(mio.num_points(r.mesh), 5);
    expect_equal(sort(r.node_perm), (1:5).');
    expect_equal(numel(r.cell_perms), 1);
    expect_equal(sort(r.cell_perms{1}), [1; 2]);
end

function t_decimate()
    f = fan();
    d = mio.decimate(f, 'target_faces', 1);
    expect_equal(mio.num_cells(d.mesh), 2);
    expect_equal(d.faces_removed, 2);
    expect_equal(d.points_removed, 1); % the unfrozen centre collapses
    assert(iscell(d.cell_maps));
    % frozen: pinning the centre vertex (id 5) prevents its own removal.
    df = mio.decimate(f, 'target_faces', 1, 'frozen', 5);
    expect_equal(df.points_removed, 0);
    % An out-of-range frozen id fails by name (0-based internally).
    expect_error(@() mio.decimate(f, 'target_faces', 1, 'frozen', 999), 'frozen node id 998');
    dv = mio.decimate(f, 'max_deviation', 1e-9);
    expect_equal(dv.deviation_met, true);
    expect_error(@() mio.decimate(f, 'max_deviation', 0.1, 'ratio', 0.5), 'max_deviation');
end

function t_error_size()
    m = fixture();
    if has_cmd('mio_estimate_error') && exist('mio.estimate_error', 'file')
        e = mio.estimate_error(m, 'temperature');
        src = e.mesh;
    else
        % The indicator mio.estimate_error would write, attached by hand.
        src = fixture();
        mio.append_cell_data(src, 'error:zz', [0.1; 0.2]);
    end
    s = mio.error_size(src, 1e-3);
    assert(any(strcmp('size', mio.point_data_names(s.mesh))));
    assert(s.num_cells + s.num_skipped >= 1);
    assert(s.min_size <= s.max_size);
    t = mio.error_size(src, 1e-3, 'output', 'h', 'min_size', 0.5, 'max_size', 2);
    assert(any(strcmp('h', mio.point_data_names(t.mesh))));
    % A target error that is not positive, crossed bounds and a missing
    % indicator are refused by name.
    expect_error(@() mio.error_size(src, 0));
    expect_error(@() mio.error_size(src, 1e-3, 'min_size', 2, 'max_size', 1));
    expect_error(@() mio.error_size(src, 1e-3, 'error_array', 'nope'));
end

function t_check_tetrahedra()
    m = fixture();
    r = mio.check_tetrahedra(m);
    assert(r.num_tetrahedra == 2);
    assert(r.num_inverted + r.num_degenerate + r.num_overlaps >= 0);
    assert(r.min_volume <= r.max_volume);
    assert(islogical(r.valid));
    assert(iscellstr(r.warnings));
    r = mio.check_tetrahedra(m, 'zero_volume_factor', 1e-6);
    assert(r.num_tetrahedra == 2);
    expect_error(@() mio.check_tetrahedra(m, 'zero_volume_factor', -1));
end

function t_conformity()
    m = fixture();
    c = mio.check_conformity(m);
    assert(c.dimension >= 2);
    assert(c.conforming);
    expect_equal(c.interface_mode, 'conforming');
    assert(any(strcmp('conformity:hanging', mio.point_data_names(c.mesh))));
    assert(iscellstr(c.warnings));
    expect_error(@() mio.check_conformity(m, 'tolerance', -1));
end

function t_merge_interpolate()
    a = fixture();
    b = fixture();
    merged = mio.merge({a, b});
    expect_equal(mio.num_points(merged), 10);
    expect_equal(mio.num_cells(merged), 4);
    assert(any(strcmp('source_mesh_id', mio.cell_data_names(merged))));
    merged2 = mio.merge({a, b}, 'weld', true, 'atol', 1e-9, 'source_tag', false);
    expect_equal(mio.num_points(merged2), 5); % the two copies weld together
    expect_error(@() mio.merge({}), 'non-empty');
    expect_error(@() mio.merge({a, b}, 'data_policy', 'union'), 'data_policy');
    merged3 = mio.merge({a, b}, 'data_policy', 'fill');
    expect_equal(mio.num_cells(merged3), 4);

    % The target already carries "temperature", so the default on_conflict
    % ('error') must refuse rather than silently overwrite.
    target = fixture();
    expect_error(@() mio.interpolate(a, target, 'arrays', 'temperature'));
    it = mio.interpolate(a, target, 'arrays', 'temperature', 'on_conflict', 'suffix');
    assert(any(strcmp('temperature_interp', mio.point_data_names(it))));

    % Conservative transfer of cell data onto an identical target.
    ct = mio.conservative_interpolate(a, target, 'arrays', {'material'}, ...
                                      'on_conflict', 'suffix');
    assert(isa(ct, 'mio.MeshHandle'));
end

function t_diff()
    a = fixture();
    b = fixture();
    rep = mio.diff(a, b);
    expect_equal(rep.verdict, 'identical');
    expect_equal(rep.points.num_exceeding, 0);
    expect_equal(size(rep.blocks.type_mismatch), [1, 1]);
    assert(~rep.blocks.type_mismatch);

    moved = fixture_points();
    moved(1, 1) = 5.0;
    mio.set_points(b, moved);
    rep2 = mio.diff(a, b);
    expect_equal(rep2.verdict, 'different');
    expect_equal(rep2.points.max_abs_error, 5.0);
    assert(~mio.meshes_equal(a, b));
    assert(mio.meshes_equal(a, b, 'atol', 10)); % inside a loose tolerance
    assert(mio.meshes_equal(a, a));
end

function t_blend()
    n = 5;
    b0 = fixture();
    b1 = fixture();
    mio.add_point_data(b0, 'u', zeros(n, 1));
    mio.add_point_data(b1, 'u', 10 * ones(n, 1));
    h = mio.blend_steps(b0, b1, 0.25);
    expect_equal(double(mio.point_data(h, 'u')), 2.5 * ones(n, 1));
    expect_equal(mio.points(h), fixture_points());
    moved = fixture_points() + 1;
    mio.set_points(b1, moved);
    hp = mio.blend_steps(b0, b1, 0.5, 'blend_points', true);
    expect_equal(mio.points(hp), fixture_points() + 0.5);
end

function t_remesh_quads()
    qb = mio.remesh_quads_backends();
    assert(iscellstr(qb));
    expect_equal(size(qb, 2), 1);
    m = cube_surface();
    expect_error(@() mio.remesh_quads(m, 'backend', 'no-such-remesher'), 'no-such-remesher');
    if any(strcmp('quadriflow', qb)) && has_cmd('mio_icosa_sphere')
        s = mio.icosa_sphere('subdivisions', 3);
        rq = mio.remesh_quads(s, 'backend', 'quadriflow', 'target_faces', 400);
        expect_equal(rq.backend, 'quadriflow');
        assert(rq.num_quads > 0);
        expect_equal(rq.num_triangles, 0);
    end
end

function t_cut()
    m = cube_surface();
    % The cube skin cut by the plane x = 0.3 -- four strips of 0.3 and the
    % face x = 0 -- or the refusal naming the build flag.
    plane = {[0.3, 0, 0, 1, 0, 0]};
    if mio.has_cutcells()
        cl = mio.cut_level_set(m, 'planes', plane);
        expect_equal(cl.backend, 'lut');
        expect_equal(cl.measure, 2.2, 1e-12);
        cq = mio.cut_quadrature(m, 'planes', plane, 'selection', 'plane0 = 0');
        expect_equal(cq.sum_weights, 4, 1e-12);
        % A plane matrix (one plane per row) is the same plane.
        cm = mio.cut_level_set(m, 'planes', [0.3, 0, 0, 1, 0, 0]);
        expect_equal(cm.measure, 2.2, 1e-12);
    else
        expect_error(@() mio.cut_level_set(m, 'planes', plane), 'MESHIOPLUSPLUS_WITH_CUTCELLS');
        expect_error(@() mio.cut_quadrature(m, 'planes', plane), 'MESHIOPLUSPLUS_WITH_CUTCELLS');
    end
    expect_error(@() mio.cut_level_set(m, 'planes', [1 2 3 4 5]), 'six numbers');
end

function t_boundary_layer()
    s = grid2d(3, 3, 'triangle');
    bl = mio.boundary_layer(s, 0.01, 'num_layers', 2);
    assert(bl.num_wedges > 0);
    assert(bl.min_scaled_jacobian > 0);
    expect_equal(bl.num_fewer_layers, 0);
    assert(bl.num_front_crossings >= 0);
    bl = mio.boundary_layer(s, 0.01, 'num_layers', 2, 'collapse', 'layers', ...
                            'normal_smooth_iterations', 2, 'front_check', false);
    assert(bl.num_wedges > 0);
    expect_equal(bl.num_front_crossings, 0);
    expect_error(@() mio.boundary_layer(s, 0.01, 'collapse', 'rows'), 'collapse');
    expect_error(@() mio.boundary_layer(s, 0.01, 'normal_smooth_iterations', 5000), ...
                 'normal_smooth_iterations');
end

function t_sweeps()
    r = grid2d(3, 2, 'quad');
    e = mio.extrude(r, 'num_layers', 4);
    expect_equal(e.num_hexahedra, 24);
    expect_equal(e.num_flipped, 0);
    assert(iscellstr(e.warnings));
    g = mio.extrude(r, 'offsets', [0.25, 1]);
    expect_equal(g.num_cells, 12);
    expect_error(@() mio.extrude(r, 'direction', [0 1 0]), 'parallel');
    v = mio.revolve(r, 'num_layers', 8);
    assert(v.full);
    expect_equal(v.num_on_axis, 3);
    expect_equal(v.num_wedges, 16);
    expect_error(@() mio.revolve(r, 'num_layers', 2), '3 steps');
    if has_cmd('mio_box') && exist('mio.box', 'file')
        b = mio.box('dims', [2 2 2], 'cell_type', 'tetra');
    else
        c = mio.convert_cells(hex_box(2), 'simplexify');
        b = c.mesh;
    end
    nb = mio.num_cells(b);
    m = mio.mirror(b);
    expect_equal(m.num_flipped_cells, nb);
    st = mio.stats(m.mesh);
    expect_equal(st.num_inverted, 0);
    k = mio.mirror(b, 'keep_original', true);
    expect_equal(k.num_points_on_plane, 9);
    expect_equal(mio.num_cells(k.mesh), 2 * nb);
    expect_error(@() mio.mirror(b, 'normal', [0 0 0]), 'normal');
end

function t_mmg()
    sq = grid2d(4, 4, 'triangle');
    if mio.has_mmg()
        r = mio.adapt(sq, 'hsiz', 0.1, 'frozen', 13);
        expect_equal(r.status, 'success');
        expect_equal(r.library, 'mmg2d');
        assert(r.num_elements_out > r.num_elements_in);
        expect_equal(size(r.quality), [3, 1]);
        xy = mio.points(sq);
        mio.add_point_data(sq, 'phi', sqrt((xy(:, 1) - 0.5).^2 + (xy(:, 2) - 0.5).^2) - 0.3);
        ls = mio.discretize_levelset(sq, 'phi', 'hmax', 0.1);
        expect_equal(ls.library, 'mmg2d');
        box = mio.mesh();
        mio.set_points(box, [0 0; 1 0; 1 1; 0 1]);
        mio.add_cell_block(box, 'line', [1 2; 2 3; 3 4; 4 1]);
        d = mio.mesh_domain(box, 'hsiz', 0.1);
        expect_equal(d.num_subdomains, 1);
        expect_equal(d.num_holes, 0);
        assert(d.num_elements_out > 100);
        expect_error(@() mio.move(sq, 'phi', 'nope'));
    else
        assert(~mio.has_mmg_elas());
        expect_error(@() mio.adapt(sq, 'hsiz', 0.1), 'MESHIOPLUSPLUS_WITH_MMG');
        expect_error(@() mio.mesh_domain(sq, 'hsiz', 0.1), 'MESHIOPLUSPLUS_WITH_MMG');
        expect_error(@() mio.discretize_levelset(sq, 'phi', 'materials', [1 1 2 3]), ...
                     'MESHIOPLUSPLUS_WITH_MMG|phi');
        expect_error(@() mio.move(sq, 'phi', {'nope'}), 'MESHIOPLUSPLUS_WITH_MMG|phi|nope');
    end
end

function t_metric()
    m = fixture();
    mio.add_point_data(m, 'hess', repmat([2 0 0 0 2 0 0 0 2], 5, 1));
    out = mio.metric_from_hessian(m, 'hess');
    assert(any(strcmp('metric', mio.point_data_names(out))));
    expect_equal(mio.num_points(out), 5);
    expect_error(@() mio.metric_from_hessian(m, 'nope'));
end

function t_tetrahedralize()
    if has_cmd('mio_icosa_sphere') && exist('mio.icosa_sphere', 'file')
        s = mio.icosa_sphere('subdivisions', 2);
        sphere = true;
    else
        s = cube_surface();
        sphere = false;
    end
    % A mesher name no program answers to is refused naming it.
    expect_error(@() mio.tetrahedralize(s, 'backend', 'no-such-mesher'), 'no-such-mesher');
    expect_error(@() mio.tetrahedralize(s, 'grading', 2), 'grading');
    if mio.has_netgen()
        t = mio.tetrahedralize(s, 'max_size', 0.3);
        assert(t.num_tetrahedra > 0);
        expect_equal(t.num_regions, 1);
        assert(t.boundary_preserved);
        expect_equal(t.backend, 'netgen');
        expect_equal(t.num_points - t.num_new_points, mio.num_points(s));
        if sphere
            expect_equal(t.num_boundary_triangles, 320);
            assert(t.volume > 3.9 && t.volume < 4 / 3 * pi);
        end
    else
        t = mio.tetrahedralize(s);
        expect_equal(t.backend, 'octree');
        expect_error(@() mio.tetrahedralize(s, 'backend', 'netgen'), 'MESHIOPLUSPLUS_WITH_NETGEN');
    end
    % The native octree backend, always there: the boundary fitted.
    t = mio.tetrahedralize(s, 'backend', 'octree');
    expect_equal(t.backend, 'octree');
    assert(~t.boundary_preserved);
    assert(t.num_tetrahedra > 0);
    if sphere
        assert(t.hausdorff_to_surface > 0);
    else
        assert(t.hausdorff_to_surface >= 0);
    end
end

function t_intersections()
    if has_cmd('mio_icosa_sphere') && exist('mio.icosa_sphere', 'file')
        a = mio.icosa_sphere('subdivisions', 3);
        b = mio.icosa_sphere('subdivisions', 3, 'center', [1 0 0]);
        wall = 2;
        tol = 0.05;
    else
        a = cube_surface();
        b = cube_surface();
        mio.set_points(b, mio.points(a) + 0.5);
        wall = 1;
        tol = 1e-9;
    end
    one = mio.find_intersections(a);
    assert(one.clean);
    assert(isempty(one.other));
    assert(iscellstr(one.warnings));
    two = mio.find_intersections(a, 'other', b);
    assert(two.num_pairs > 0);
    expect_equal(size(two.pairs), [two.num_pairs, 2]);
    assert(min(two.pairs(:)) >= 1);
    assert(isa(two.other, 'mio.MeshHandle'));
    t = mio.thickness(a);
    expect_equal(t.direction, 'inward');
    expect_equal(t.num_missed, 0);
    assert(abs(t.min - wall) < tol);
    expect_equal(size(t.min_point), [1, 3]);
    assert(t.min_cell >= 1);
    expect_error(@() mio.thickness(a, 'direction', 'sideways'), 'direction');
    ts = mio.thickness(a, 'method', 'cone', 'two_sided', true);
    expect_equal(ts.method, 'cone');
    assert(ts.p2 <= ts.p95);
    expect_error(@() mio.thickness(a, 'method', 'laser'), 'method');
    ch = mio.slice(a, [0 0 0.1], [0 0 1], 'chain', true);
    assert(mio.num_points(ch) > 0);
    ss = mio.sections(a, [0 0 1], [-0.5 0.5]);
    assert(mio.num_points(ss) > 0);
    fc = mio.fit_curves(ch);
    if wall == 2
        % The sphere's section at z = 0.1 is one circle.
        expect_equal(fc.num_circles, 1);
        expect_equal(fc.pieces(1).type, 'circle');
    else
        % The cube's section is a square: four lines, never a circle through
        % its corners.
        expect_equal(fc.num_circles, 0);
        expect_equal(fc.num_lines, 4);
        expect_equal(fc.pieces(1).type, 'line');
    end
    sil = mio.silhouette(a, 'clearance', 0.05);
    expect_equal(sil.num_loops, 1);
    expect_equal(sil.num_holes, 0);
    assert(sil.area > 0);
    expect_error(@() mio.silhouette(a, 'arc_segments', 2), 'arc_segments');
end

function t_capabilities()
    caps = {mio.has_netgen(), mio.has_cutcells(), mio.has_mmg(), mio.has_mmg_elas()};
    for k = 1:numel(caps)
        assert(islogical(caps{k}) && isscalar(caps{k}));
    end
    tb = mio.tetrahedralize_backends();
    assert(iscellstr(tb));
    assert(any(strcmp('octree', tb)));
    expect_equal(size(tb, 2), 1);
end

function t_decimate_volume()
    tets = mio.box('dims', [3 3 3], 'cell_type', 'tetra');
    n_in = mio.num_cells(tets);
    d = mio.decimate_volume(tets, 'ratio', 0.5);
    assert(d.tets_removed > 0);
    expect_equal(mio.num_cells(d.mesh), n_in - d.tets_removed);
    expect_equal(mio.cell_block_types(d.mesh), {'tetra'});
    % 1-based maps, 0 for an entity that did not survive.
    expect_equal(size(d.point_map), [mio.num_points(tets), 1]);
    cm = d.cell_maps{1};
    expect_equal(size(cm), [n_in, 1]);
    assert(all(cm >= 0));
    expect_equal(sum(cm > 0), mio.num_cells(d.mesh));
    % Exactly one criterion, and frozen ids are range-checked.
    expect_error(@() mio.decimate_volume(tets, 'ratio', 0.5, 'target_cells', 10));
    expect_error(@() mio.decimate_volume(tets, 'ratio', 0.5, 'frozen', 99999), 'frozen');
end

function t_deviation()
    if ~(has_cmd('mio_icosa_sphere') && exist('mio.icosa_sphere', 'file'))
        return;
    end
    ref = mio.icosa_sphere('subdivisions', 3);
    scan = mio.icosa_sphere('radius', 1.01, 'subdivisions', 3);
    d = mio.deviation(scan, 'reference', ref, 'tolerance', 0.02);
    assert(abs(d.mean - 0.01) < 1e-9);
    expect_equal(d.share_within, 1);
    expect_equal(numel(d.histogram_counts), 20);
    expect_equal(sum(d.histogram_counts), d.num_measured - d.num_outside_range);
    p = mio.deviation(scan, 'type', 'sphere', 'params', [0 0 0 0 0 0 1 0 0 0]);
    assert(abs(p.mean - 0.01) < 1e-12);
    expect_error(@() mio.deviation(scan, 'type', 'blob'), 'type');
    g = mio.deviation_regions(d.mesh, 0.005);
    expect_equal(g.num_excess, 1);
    expect_equal(g.regions(1).name, 'deviation:excess_0');
    e = mio.exaggerate(d.mesh, 'scale', 20);
    assert(e.num_moved > 0);
    c = mio.colorize(d.mesh, 'deviation:distance', 'steps', 4);
    expect_equal(c.location, 'point');
end
