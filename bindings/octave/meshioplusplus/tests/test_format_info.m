function tests = test_format_info()
%TEST_FORMAT_INFO  Format side channels, the isogeometric layer and CAD
%   meshing (format_info.c). Ports test-gmsh-periodic.R, test-mfem-nurbs.R,
%   test-cad-json-brep.R and the side-channel / mesh_cad_surface tests of
%   test-meshioplusplus.R.
    tests = {
        'Gmsh periodic links use the format side channel', @t_gmsh_periodic
        'an MFEM NURBS mesh keeps its patches in the side channel', @t_mfem_nurbs
        'a Kratos CAD JSON file keeps its B-rep in the side channel', @t_cad_json_brep
        'the isogeometric layer refines a trimmed face and tests inside', @t_iga_trimmed
        'an mdpa side channel survives read_with_info / write_with_info', @t_mdpa
        'mesh_cad_surface meshes a trimmed face and reports it', @t_mesh_cad_surface
        'format info handles display, release and reject misuse', @t_handles
    };
end

% --- fixtures ---------------------------------------------------------------

function write_lines(path, lines)
    fid = fopen(path, 'w');
    assert(fid > 0, 'could not open %s', path);
    for k = 1:numel(lines)
        fprintf(fid, '%s\n', lines{k});
    end
    fclose(fid);
end

function line = first_line(path)
    fid = fopen(path, 'r');
    line = fgetl(fid);
    fclose(fid);
end

function path = mfem_fixture()
    path = [tempname() '.mesh'];
    w = '0.70710678118654757';
    write_lines(path, {'MFEM NURBS mesh v1.0', 'dimension', '2', 'elements', '1', ...
        '5 3 0 1 2 3', 'boundary', '0', 'edges', '4', '0 0 1', '0 3 2', '1 0 3', '1 1 2', ...
        'vertices', '4', 'knotvectors', '2', '2 3 0 0 0 1 1 1', '1 2 0 0 1 1', ...
        'weights', '1', '1', '1', '1', w, w, ...
        'FiniteElementSpace', 'FiniteElementCollection: NURBS2', 'VDim: 2', 'Ordering: 1', ...
        '1 0', '0 1', '0 2', '2 0', '1 1', '2 2'});
end

% The trimmed triangle face of test-cad-json-brep.R and the mesh_cad_surface
% test: a 2 x 2 bilinear patch trimmed to the triangle (0,0), (1.5,0), (0,1.5).
function path = triangle_cad_fixture()
    path = [tempname() '.cad.json'];
    write_lines(path, {[ ...
        '{"breps": [{"brep_id": 1, "faces": [{"brep_id": 2, "surface": {"is_trimmed": true, ' ...
        '"degrees": [1, 1], "knot_vectors": [[0, 2], [0, 2]], "control_points": ' ...
        '[[1, [0, 0, 0, 1]], [2, [2, 0, 0, 1]], [3, [0, 2, 0, 1]], [4, [2, 2, 0, 1]]]}, ' ...
        '"boundary_loops": [{"loop_type": "outer", "trimming_curves": [' ...
        '{"trim_index": 0, "curve_direction": true, "parameter_curve": {"degree": 1, ' ...
        '"knot_vector": [0, 0, 1, 1], "active_range": [0, 1], "control_points": [[[0, 0, 0, 1]], [[1.5, 0, 0, 1]]]}}, ' ...
        '{"trim_index": 1, "curve_direction": true, "parameter_curve": {"degree": 1, ' ...
        '"knot_vector": [0, 0, 1, 1], "active_range": [0, 1], "control_points": [[[1.5, 0, 0, 1]], [[0, 1.5, 0, 1]]]}}, ' ...
        '{"trim_index": 2, "curve_direction": true, "parameter_curve": {"degree": 1, ' ...
        '"knot_vector": [0, 0, 1, 1], "active_range": [0, 1], "control_points": [[[0, 1.5, 0, 1]], [[0, 0, 0, 1]]]}}]}]}], ' ...
        '"edges": []}]}']});
end

% --- tests ------------------------------------------------------------------

function t_gmsh_periodic()
    path = [tempname() '.msh'];
    output = [tempname() '.msh'];
    cleanup = onCleanup(@() delete_files({path, output}));
    write_lines(path, {'$MeshFormat', '4.0 0 8', '$EndMeshFormat', ...
        '$Nodes', '1 3', '22 2 0 3', '30 1 0 0', '10 0 0 0', '20 0 1 0', ...
        '$EndNodes', '$Elements', '1 1', '22 2 2 1', '90 10 30 20', '$EndElements', ...
        '$Periodic', '1', '1 12 33', '3', '30 10', '20 10', '30 10', '$EndPeriodic'});
    expect_error(@() mio.read(path, 'format', 'gmsh'), 'Periodic');
    r = mio.read_with_info(path, 'format', 'gmsh');
    copied = mio.gmsh_info(r.info);
    expect_equal(copied.bounding_entities, cell(0, 1));
    expect_equal(numel(copied.periodic), 1);
    expect_equal(copied.periodic(1).entity, [1; 12; 33]);
    expect_equal(copied.periodic(1).affine, zeros(0, 1));
    % R's (2, N) node_pairs are (N x 2) here: one (slave, master) per row.
    expect_equal(copied.periodic(1).node_pairs, [1, 2; 3, 2; 1, 2]);
    formats = {'gmsh', 'gmsh22'};
    for k = 1:numel(formats)
        mio.write_with_info(r.mesh, r.info, output, 'format', formats{k});
        back = mio.read_with_info(output, 'format', 'gmsh');
        g = mio.gmsh_info(back.info);
        expect_equal(g.periodic, copied.periodic);
        mio.release(back.mesh);
        mio.format_info_release(back.info);
    end
    expect_error(@() mio.mdpa_info(r.info));
    mio.release(r.mesh);
    mio.format_info_release(r.info);
    expect_equal(copied.periodic(1).node_pairs(2, 2), 2);
    expect_error(@() mio.gmsh_info(r.info), 'released');
end

function t_mfem_nurbs()
    path = mfem_fixture();
    output = [tempname() '.mesh'];
    refined_path = [tempname() '.mesh'];
    cleanup = onCleanup(@() delete_files({path, output, refined_path}));
    r = mio.read_with_info(path, 'format', 'mfem');
    expect_equal(mio.format_info_format(r.info), 'mfem');
    copied = mio.nurbs_info(r.info);
    expect_equal(numel(copied.patches), 1);
    p = copied.patches(1);
    expect_equal(p.degrees, [2; 1]);
    expect_equal(p.knots{2}, [0; 0; 1; 1]);
    expect_equal(size(p.control_points), [6, 2]);
    expect_equal(p.control_points(2, :), [1, 1]);
    expect_equal(p.weights(2), sqrt(0.5), 1e-12);
    expect_equal(copied.field_names, cell(0, 1));
    expect_equal(copied.topology.num_vertices, 4);
    expect_equal(copied.topology.patch_vertices(1, :), [1, 2, 3, 4]);
    expect_equal(copied.topology.edges(1, :), [0, 1, 2]);
    expect_equal(copied.cell_patch, 1);
    mio.write_with_info(r.mesh, r.info, output, 'format', 'mfem');
    expect_equal(first_line(output), 'MFEM NURBS mesh v1.0');
    back = mio.read_with_info(output, 'format', 'mfem');
    nb = mio.nurbs_info(back.info);
    expect_patches_equal(nb.patches, copied.patches);
    mio.release(back.mesh);
    mio.format_info_release(back.info);
    % the isogeometric layer: refine, write back, evaluate
    fine = mio.iga_refine(r.info, 'subdivisions', 3, 'elevate', 1);
    expect_equal(mio.format_info_format(fine.info), 'mfem');
    nf = mio.nurbs_info(fine.info);
    expect_equal(nf.patches(1).degrees, [3; 2]);
    mio.write_with_info(fine.mesh, fine.info, refined_path, 'format', 'mfem');
    expect_equal(first_line(refined_path), 'MFEM NURBS mesh v1.0');
    % R's (d, N) matrix(c(0.3, 0.7, 0.9, 0.2), nrow = 2) is one point per row.
    uv = [0.3, 0.7; 0.9, 0.2];
    a = mio.iga_evaluate(r.info, 1, uv, 'derivatives', 1);
    b = mio.iga_evaluate(fine.info, 1, uv);
    expect_equal(size(a.points), [2, 2]);
    assert(max(abs(a.points(:) - b.points(:))) < 1e-13);
    expect_equal(size(a.derivatives), [2, 2, 2]);
    assert(isempty(a.normals));
    assert(isempty(a.inside));
    expect_equal(numel(a.fields), 0);
    assert(isempty(b.derivatives));
    expect_error(@() mio.iga_evaluate(r.info, 1, [2, 0.5]), 'knot range');
    expect_error(@() mio.iga_evaluate(r.info, 1, [0.5, 0.5, 0.5]), 'params');
    sampled = mio.iga_sample(r.info);
    expect_equal(mio.format_info_format(sampled.info), 'mfem');
    assert(mio.num_points(sampled.mesh) > 0);
    mio.release(sampled.mesh);
    mio.format_info_release(sampled.info);
    mio.release(fine.mesh);
    mio.format_info_release(fine.info);
    expect_error(@() mio.gmsh_info(r.info));
    mio.release(r.mesh);
    mio.format_info_release(r.info);
end

function t_cad_json_brep()
    path = [tempname() '.cad.json'];
    output = [tempname() '.cad.json'];
    cleanup = onCleanup(@() delete_files({path, output}));
    write_lines(path, { ...
        '{"breps": [{"brep_id": 1, "faces": [{"brep_id": 2, "surface": {', ...
        '"is_trimmed": true, "degrees": [1, 1], "knot_vectors": [[0, 1], [0, 1]],', ...
        '"control_points": [[1, [0, 0, 0, 1]], [2, [2, 0, 0, 1]], [3, [0, 2, 0, 1]], [4, [2, 2, 0, 1]]]},', ...
        '"boundary_loops": [{"loop_type": "outer", "trimming_curves": [', ...
        '{"trim_index": 0, "curve_direction": true, "parameter_curve": {"degree": 1,', ...
        '"knot_vector": [0, 0, 1, 1], "active_range": [0, 1],', ...
        '"control_points": [[[0, 0, 0, 1]], [[1, 0, 0, 1]]]}}]}]}],', ...
        '"edges": [{"brep_id": 3, "topology": [{"brep_id": 2, "trim_index": 0}]}]}]}'});
    r = mio.read_with_info(path);
    expect_equal(mio.format_info_format(r.info), 'cad_json');
    b = mio.brep_info(r.info);
    expect_equal(numel(b.faces), 1);
    f = b.faces(1);
    expect_equal(f.id, 2);
    expect_equal(f.patch, 1);
    assert(f.trimmed);
    assert(f.reduced_knots);
    expect_equal(f.control_point_ids, [1; 2; 3; 4]);
    expect_equal(b.patches(1).knots{1}, [0; 0; 1; 1]);
    expect_equal(size(b.patches(1).control_points), [4, 3]);
    expect_equal(f.loops(1).type, 'outer');
    expect_equal(f.loops(1).trims(1).trim_index, 0);
    assert(isempty(b.edges(1).curve));
    expect_equal(b.edges(1).topology(1).relative_direction, -1);
    expect_equal(b.breps(1).edges, 1);
    assert(ischar(b.extra));
    mio.write_with_info(r.mesh, r.info, output);
    back = mio.read_with_info(output);
    bb = mio.brep_info(back.info);
    expect_patches_equal(bb.patches, b.patches);
    mio.release(back.mesh);
    mio.format_info_release(back.info);
    expect_error(@() mio.nurbs_info(r.info));
    mio.release(r.mesh);
    mio.format_info_release(r.info);
end

function t_iga_trimmed()
    path = triangle_cad_fixture();
    cleanup = onCleanup(@() delete_files({path}));
    r = mio.read_with_info(path);
    fine = mio.iga_refine(r.info, 'subdivisions', [2, 3], 'elevate', [0, 1]);
    bf = mio.brep_info(fine.info);
    expect_equal(bf.patches(1).degrees, [1; 2]);
    uv = [0.2, 0.2; 1.8, 1.8];
    ev = mio.iga_evaluate(fine.info, 1, uv, 'derivatives', 1);
    expect_equal(ev.inside, [true; false]);
    expect_equal(abs(ev.normals(1, 3)), 1, 1e-13);
    % R's (3, 2, 2) = (sdim, d, N) is (N x d x sdim) here.
    expect_equal(size(ev.derivatives), [2, 2, 3]);
    sampled = mio.iga_sample(r.info, 'trim_samples', 8, 'edges', false);
    assert(mio.num_points(sampled.mesh) > 0);
    mio.release(sampled.mesh);
    mio.format_info_release(sampled.info);
    mio.release(fine.mesh);
    mio.format_info_release(fine.info);
    mio.release(r.mesh);
    mio.format_info_release(r.info);
end

function t_mdpa()
    dir = tempname();
    mkdir(dir);
    cleanup = onCleanup(@() rmdir(dir, 's'));
    path = fullfile(dir, 'side.mdpa');
    write_lines(path, { ...
        'Begin ModelPartData', '    SOLVER_TYPE static', 'End ModelPartData', ...
        'Begin Table 1 TIME VALUE', '    0.0 1.0', '    2.0 3.0', 'End Table', ...
        'Begin Nodes', '1 0 0 0', '2 1 0 0', '3 0 1 0', '4 0 0 1', 'End Nodes', ...
        'Begin Elements Element3D4N', '1 0 1 2 3 4', 'End Elements', ...
        'Begin Geometries Triangle3D3', '7 2 3 4', 'End Geometries', ...
        'Begin SubModelPart Part', '    Begin SubModelPart Inner', ...
        '        Begin SubModelPartGeometries', '            7', '        End SubModelPartGeometries', ...
        '        Begin SubModelPartConstraints', '            1', '        End SubModelPartConstraints', ...
        '    End SubModelPart', 'End SubModelPart', ...
        'Begin Mesh 5', '    Begin MeshNodes', '        4', '    End MeshNodes', 'End Mesh', ...
        'Begin Constraints LinearMasterSlaveConstraint', ...
        '    1 1 DISPLACEMENT_X 2 DISPLACEMENT_X 1.0 0.0', 'End Constraints'});

    expect_error(@() mio.read(path), 'mio_read_with_info');
    r = mio.read_with_info(path);
    assert(isa(r.info, 'mio.FormatInfo'));
    expect_equal(mio.format_info_format(r.info), 'mdpa');
    d = mio.mdpa_info(r.info);
    expect_equal(d.model_part_data.SOLVER_TYPE, 'static');
    expect_equal(d.tables(1).header, '1 TIME VALUE');
    expect_equal(d.tables(1).values, [0, 1; 2, 3]);
    % R's (nodes_per_geometry, num_geometries) is one geometry per row here.
    expect_equal(d.geometries(1).connectivity, [2, 3, 4]);
    expect_equal(d.geometries(1).ids, 7);
    expect_equal(d.submodelparts(1).name, 'Part/Inner');
    expect_equal(d.submodelparts(1).geometry_ids, 7);
    expect_equal(d.submodelparts(1).constraint_ids, 1);
    expect_equal(d.mesh_blocks(1).id, 5);
    expect_equal(d.mesh_blocks(1).nodes, 4);
    expect_equal(d.raw_blocks(1).terminator, 'End Constraints');
    assert(iscell(d.skipped));
    assert(isstruct(d.properties) && isstruct(d.entity_names));

    out = fullfile(dir, 'out.mdpa');
    mio.write_with_info(r.mesh, r.info, out);
    r2 = mio.read_with_info(out);
    d2 = mio.mdpa_info(r2.info);
    expect_equal(d2.raw_blocks, d.raw_blocks);
    expect_equal(d2.submodelparts, d.submodelparts);
    expect_error(@() mio.write_with_info(r.mesh, r.info, fullfile(dir, 'out.vtu')));
    mio.release(r2.mesh);
    mio.format_info_release(r2.info);
    mio.format_info_release(r2.info); % idempotent
    assert(~mio.format_info_is_open(r2.info));

    vtu = fullfile(dir, 'plain.vtu');
    mio.write(r.mesh, vtu);
    r3 = mio.read_with_info(vtu);
    assert(isempty(r3.info));
    % info = [] makes write_with_info exactly mio.write.
    mio.write_with_info(r3.mesh, [], fullfile(dir, 'plain2.vtu'));
    assert(exist(fullfile(dir, 'plain2.vtu'), 'file') == 2);
    mio.release(r3.mesh);
    mio.release(r.mesh);
    mio.format_info_release(r.info);
end

function t_mesh_cad_surface()
    path = triangle_cad_fixture();
    cleanup = onCleanup(@() delete_files({path}));
    r = mio.read_with_info(path);
    out = mio.mesh_cad_surface(r.info, 0.1);
    expect_equal(out.num_faces, 1);
    expect_equal(out.num_failed, 0);
    assert(out.num_triangles > 100);
    expect_equal(out.num_points, mio.num_points(out.mesh));
    assert(out.min_angle_deg > 10);
    expect_equal(out.face_status, {'ok'});
    expect_equal(out.face_reason, {''});
    expect_equal(numel(out.warnings), 0);
    assert(iscell(out.warnings));
    expect_error(@() mio.mesh_cad_surface(r.info, 0));
    expect_equal(out.min_input_angle_deg, 180);
    expect_equal(numel(out.corners.face), 0);
    % v17.15.0: refined to 25 degrees and smoothed; the two 45-degree corners
    % are reported.
    q = mio.mesh_cad_surface(r.info, 0.1, 'min_angle', 25, 'smooth_iterations', 2, ...
        'smooth_method', 'odt');
    assert(q.min_angle_deg >= 25 - 1e-9);
    assert(q.num_smoothed > 0);
    expect_equal(q.min_input_angle_deg, 45, 1e-9);
    expect_equal(numel(q.angle_histogram), 18);
    expect_equal(sum(q.angle_histogram), 3 * q.num_triangles);
    expect_equal(numel(q.corners.face), 2);
    assert(all(q.corners.face == 1 & ~q.corners.pole));
    expect_equal(q.corners.angle_deg, [45; 45], 1e-9);
    mio.release(q.mesh);
    expect_error(@() mio.mesh_cad_surface(r.info, 0.1, 'smooth_method', 'none'));
    expect_error(@() mio.mesh_cad_surface(r.info, 0.1, 'min_angle', 60));
    % v17.16.0: curved triangles, their nodes on the geometry.
    for o = 2:3
        c2 = mio.mesh_cad_surface(r.info, 0.1, 'order', o);
        expect_equal(c2.num_triangles, out.num_triangles);
        assert(c2.num_points > out.num_points);
        expect_equal(c2.num_curved_straightened, 0);
        mio.release(c2.mesh);
    end
    expect_error(@() mio.mesh_cad_surface(r.info, 0.1, 'order', 4));
    % v17.18.0: the size field and pole rings.
    s = mio.mesh_cad_surface(r.info, 0.1, 'grading', 1.3, 'curvature_deviation', 0.01, ...
        'cells_across_gaps', 2, 'pole_rings', 2);
    expect_equal(s.num_failed, 0);
    mio.release(s.mesh);
    expect_error(@() mio.mesh_cad_surface(r.info, 0.1, 'grading', 0.5));
    % v17.19.0: a size given to face 1 (1-based) refines it. R's
    % c("1" = 0.025) is the [index, size] row [1, 0.025].
    g = mio.mesh_cad_surface(r.info, 0.1, 'face_sizes', [1, 0.025]);
    expect_equal(g.num_failed, 0);
    assert(g.num_triangles > 400);
    assert(iscell(g.face_size_limits));
    expect_equal(numel(g.virtual_edges), 0);
    mio.release(g.mesh);
    expect_error(@() mio.mesh_cad_surface(r.info, 0.1, 'face_sizes', [7, 0.025]));
    expect_error(@() mio.mesh_cad_surface(r.info, 0.1, 'face_sizes', [1; 0.025]), 'n x 2');
    % v17.20.0: a planar face is never flattened, and nothing is a sliver.
    x = mio.mesh_cad_surface(r.info, 0.1, 'reparam_condition_number', 4, ...
        'model_tolerances', true, 'virtual_faces', true);
    expect_equal(x.num_failed, 0);
    expect_equal(numel(x.reparametrised.face), 0);
    expect_equal(numel(x.merged_faces.sliver), 0);
    expect_equal(numel(x.tolerance_uses.edge), 0);
    expect_equal(numel(x.lopped_corners.face), 0);
    mio.release(x.mesh);
    % Corner lopping at 20 degrees: the two 45-degree corners, under the
    % straight-segment bound.
    l = mio.mesh_cad_surface(r.info, 0.1, 'min_angle', 20, 'lop_corners', true);
    expect_equal(numel(l.lopped_corners.face), 2);
    expect_equal(l.lopped_corners.angle_deg, [45; 45], 1e-9);
    assert(all(l.lopped_corners.guaranteed));
    assert(all(l.lopped_corners.stated_bound_deg > 20.7));
    assert(iscellstr(l.lopped_corners.reason));
    mio.release(l.mesh);
    expect_error(@() mio.mesh_cad_surface(r.info, 0.1, 'reparam_condition_number', 0.5));
    % A metric error without a metric background is refused (roadmap 6.10.2).
    expect_error(@() mio.mesh_cad_surface(r.info, 0.1, 'metric_error', 1));
    % Quads (roadmap 6.8) and untangling (6.6.3).
    q = mio.mesh_cad_surface(r.info, 0.1, 'min_angle', 25, 'quads', 'dominant');
    assert(q.num_quads > 0);
    assert(q.min_quad_scaled_jacobian > 0);
    expect_equal(q.min_metric_angle, 0);
    mio.release(q.mesh);
    g = mio.mesh_cad_surface(r.info, 0.1, 'min_angle', 25, 'quads', 'dominant', ...
        'quad_pairs', 'greedy');
    assert(g.num_quads > 0);
    mio.release(g.mesh);
    expect_error(@() mio.mesh_cad_surface(r.info, 0.1, 'untangle', true));
    u = mio.mesh_cad_surface(r.info, 0.1, 'order', 2, 'untangle', true);
    expect_equal(u.num_untangled, 0);
    mio.release(u.mesh);
    % Projection back onto the model (roadmap 6.10.1) moves nothing.
    p = mio.mesh_cad_surface(r.info, 0.1);
    b = mio.project_to_cad(p.mesh, r.info);
    expect_equal(b.num_projected, mio.num_points(p.mesh));
    assert(b.max_move < 1e-12);
    assert(b.num_on_curves >= 0);
    assert(b.num_on_curves <= b.num_on_edges);
    expect_equal(numel(b.outside_trim), 0);
    mio.release(b.mesh);
    b = mio.project_to_cad(p.mesh, r.info, 'cad:face', 'cad:edge', 'max_distance', 1, ...
        'trim_samples', 8);
    expect_equal(b.num_projected, mio.num_points(p.mesh));
    mio.release(b.mesh);
    mio.release(p.mesh);
    mio.release(out.mesh);
    mio.release(r.mesh);
    mio.format_info_release(r.info);
end

function t_handles()
    path = mfem_fixture();
    cleanup = onCleanup(@() delete_files({path}));
    r = mio.read_with_info(path, 'format', 'mfem');
    assert(mio.format_info_is_open(r.info));
    s = evalc('disp(r.info)');
    assert(~isempty(strfind(s, '<mio.FormatInfo: mfem>')), s);
    % Wrong kinds are rejected, never dereferenced.
    expect_error(@() mio.format_info_format(r.mesh), 'mio.FormatInfo');
    expect_error(@() mio.format_info_release(r.mesh), 'mio.FormatInfo');
    expect_error(@() mio.nurbs_info([]));
    assert(~mio.format_info_is_open(r.mesh));
    assert(~mio.format_info_is_open(42));
    mio.format_info_release(r.info);
    mio.format_info_release(r.info);
    assert(~mio.format_info_is_open(r.info));
    s = evalc('disp(r.info)');
    assert(~isempty(strfind(s, '<mio.FormatInfo: released>')), s);
    expect_error(@() mio.format_info_format(r.info), 'released');
    mio.release(r.mesh);
end

% --- helpers ----------------------------------------------------------------

% R compares patch lists with expect_equal's tolerance; the values pass
% through a text file, so compare numbers within a relative 1e-12.
function expect_patches_equal(a, b)
    expect_equal(size(a), size(b));
    expect_equal(fieldnames(a), fieldnames(b));
    f = fieldnames(a);
    for k = 1:numel(a)
        for j = 1:numel(f)
            x = a(k).(f{j});
            y = b(k).(f{j});
            if iscell(x)
                expect_equal(size(x), size(y));
                for c = 1:numel(x)
                    expect_equal(x{c}, y{c});
                end
            else
                expect_equal(x, y);
            end
        end
    end
end

function delete_files(paths)
    for k = 1:numel(paths)
        if exist(paths{k}, 'file') == 2
            delete(paths{k});
        end
    end
end
