function tests = test_ops_data()
%TEST_OPS_DATA  Data operations, the settings pipeline, grids, generators,
%   voxelize / octree / signed distance, remesh_volume, provenance and select
%   (ops_data.c).
    tests = {
        'region-backed sets/data conversions preserve labels and geometry', @t_sets
        'data operations never touch geometry', @t_data_ops
        'data operation options reach the C API', @t_data_ops_options
        'tensor invariants: von Mises / principal / hydrostatic / deviatoric', @t_tensor
        'field integration totals/means over cells', @t_integrate
        'the settings pipeline runs (or fails naming the flag)', @t_pipeline
        'pipeline Version 2 reaches terminal partition through report entry points', @t_pipeline_v2
        'sequence pipeline reports run a plain document', @t_sequence_report
        'grid builds a lattice from nothing', @t_grid
        'the parametric generators build meshes and refuse bad requests', @t_generators
        'every generator runs with its defaults', @t_generator_smoke
        'reference cells take a list or a string', @t_reference_cells
        'voxelize builds a grid around a surface', @t_voxelize
        'octree_mesh meshes the inside of a skin', @t_octree_mesh
        'signed distance matches the cube''s closed form', @t_distance
        'distance at cell centres', @t_distance_center
        'compute_sdf generates the grid and the field in one call', @t_compute_sdf
        'remesh_volume retetrahedralizes by isosurface stuffing', @t_remesh_volume
        'select writes a named region', @t_select
        'provenance scopes open, record and close', @t_provenance
        'errors carry the C library message', @t_errors
    };
end

% --- data operations -------------------------------------------------------

function t_sets()
    m = fixture();
    mio.add_region(m, 'a', 'point', [1 2]);
    mio.add_region(m, 'b', 'point', [2 3]);
    labels = mio.sets_to_data(m, 'point', 'order', {'b', 'a'});
    expect_equal(double(mio.point_data(labels, 'b-a')), [1; 1; 0; -1; -1]);
    expect_equal(numel(mio.regions(labels)), 0);
    expect_equal(mio.points(labels), mio.points(m));
    restored = mio.data_to_sets(labels, 'point', 'b-a');
    expect_equal(numel(mio.regions(restored)), 3);
    named = mio.sets_to_data(m, 'point', 'data_name', 'label', 'join_char', '+');
    assert(any(strcmp('label', mio.point_data_names(named))));
    expect_error(@() mio.sets_to_data(m, 'field'));
    expect_error(@() mio.data_to_sets(m, 'point', 'temperature'));
end

function t_data_ops()
    m = fixture();
    before_points = mio.points(m);
    before_conn = mio.connectivity(m, 1);

    d = mio.data_drop(m, 'point', 'temperature');
    expect_equal(mio.point_data_names(d), {'displacement'});
    expect_equal(mio.points(d), before_points); % geometry bit-identical
    expect_equal(mio.connectivity(d, 1), before_conn);

    k = mio.data_keep(m, 'point', 'temperature');
    expect_equal(mio.point_data_names(k), {'temperature'});

    r = mio.data_rename(m, 'point', 'temperature', 'T');
    assert(any(strcmp('T', mio.point_data_names(r))));
    assert(~any(strcmp('temperature', mio.point_data_names(r))));

    p2c = mio.data_point_to_cell(m, 'temperature');
    assert(any(strcmp('temperature', mio.cell_data_names(p2c))));

    c2p = mio.data_cell_to_point(m, 'material');
    assert(any(strcmp('material', mio.point_data_names(c2p))));

    calc = mio.data_calc(m, 'temperature * 2.0', 'doubled');
    expect_equal(mio.point_data(calc, 'doubled'), [2; 4; 6; 8; 10]);

    cond = mio.data_condition(m, 'point', 'names', 'temperature', 'cond_mode', 'normalize');
    t = mio.point_data(cond, 'temperature');
    expect_equal(min(t), 0);
    expect_equal(max(t), 1);

    info = mio.data_info(m);
    expect_equal(size(info), [3, 1]);
    temp = info(strcmp({info.name}, 'temperature'));
    expect_equal(temp.location, 'point');
    expect_equal(temp.dtype, 'float64');
    expect_equal(temp.num_entries, 5);
    expect_equal(temp.num_components, 1);
    expect_equal(temp.min, 1);
    expect_equal(temp.max, 5);
    expect_equal(temp.mean, 3);
    expect_equal(temp.num_nan, 0);
    expect_equal(temp.components, [1, 5, 3]);
    disp_info = info(strcmp({info.name}, 'displacement'));
    expect_equal(size(disp_info.components), [3, 3]);
    expect_equal(disp_info.components(:, 1), [10; 11; 12]); % per-component min
end

function t_data_ops_options()
    m = fixture();
    % [] and {} both mean "every array" for the averaging entry points.
    all1 = mio.data_point_to_cell(m);
    all2 = mio.data_point_to_cell(m, {}, '_c');
    assert(all(ismember({'displacement'; 'temperature'}, mio.cell_data_names(all1))));
    assert(any(strcmp('temperature_c', mio.cell_data_names(all2))));
    w = mio.data_cell_to_point(m, 'names', 'material', 'weight', 'measure', 'suffix', '_p');
    assert(any(strcmp('material_p', mio.point_data_names(w))));
    expect_error(@() mio.data_cell_to_point(m, 'weight', 'heavy'), 'weight');
    clamp = mio.data_condition(m, 'point', 'names', {'temperature'}, 'lo', 2, 'hi', 4, ...
        'suffix', '_c');
    expect_equal(mio.point_data(clamp, 'temperature_c'), [2; 2; 3; 4; 4]);
    mag = mio.data_condition(m, 'point', 'names', 'displacement', 'cond_mode', 'standardize', ...
        'scope', 'magnitude', 'nan_policy', 'replace', 'nan_replacement', 0);
    expect_equal(size(mio.point_data(mag, 'displacement')), [5, 3]);
    expect_error(@() mio.data_condition(m, 'point', 'cond_mode', 'squash'), 'cond_mode');
    expect_error(@() mio.data_drop(m, 'volume', 'x'), 'location');
    % ignore_missing skips names that do not exist.
    expect_error(@() mio.data_drop(m, 'point', 'nope'));
    kept = mio.data_drop(m, 'point', {'nope', 'temperature'}, 'ignore_missing', true);
    expect_equal(mio.point_data_names(kept), {'displacement'});
    over = mio.data_calc(m, 'temperature + 1', 'temperature', 'overwrite', true);
    expect_equal(mio.point_data(over, 'temperature'), [2; 3; 4; 5; 6]);
    expect_error(@() mio.data_calc(m, 'temperature + 1', 'temperature'));
    cell_calc = mio.data_calc(m, 'material * 2', 'm2', 'location', 'cell');
    expect_equal(mio.cell_data(cell_calc, 'm2'), [14; 18]);
end

function t_tensor()
    m = fixture();
    % xx yy zz xy yz zx, one row per point (5 points).
    stress = repmat([1, 2, 3, 0.5, 0.6, 0.7], 5, 1);
    mio.add_point_data(m, 'stress', stress);

    inv = mio.tensor_invariants(m, 'point', 'names', 'stress', ...
        'outputs', {'mises', 'hydrostatic'});
    mises = mio.point_data(inv, 'stress_mises');
    expect_equal(mises(1), ...
        sqrt(0.5 * ((1 - 2)^2 + (2 - 3)^2 + (3 - 1)^2 + 6 * (0.5^2 + 0.6^2 + 0.7^2))));
    hyd = mio.point_data(inv, 'stress_hydrostatic');
    expect_equal(hyd(1), 2);
    assert(~any(strcmp('stress_principal', mio.point_data_names(inv))));

    all4 = mio.tensor_invariants(m, 'point', 'names', 'stress', 'prefix', 'p_');
    names = mio.point_data_names(all4);
    assert(any(strcmp('p_stress_principal', names)));
    assert(any(strcmp('p_stress_deviatoric', names)));
    expect_error(@() mio.tensor_invariants(m, 'point', 'outputs', 'trace'), 'outputs');
    expect_error(@() mio.tensor_invariants(m, 'field'), 'field_data');
end

function t_integrate()
    m = fixture();
    mio.add_region(m, 'solid', 'cell', 1, 'dim', 3, 'tag', 17);

    report = mio.data_integrate(m, 'material');
    expect_equal(numel(report), 1);
    arr = report(1);
    expect_equal(arr.name, 'material');
    expect_equal(arr.num_components, 1);
    expect_equal(arr.domain.num_cells, 2);
    expect_equal(arr.domain.num_skipped, 0);
    dcomp = arr.domain.components; % columns: total, mean, domain_measure, num_nan
    expect_equal(size(dcomp), [1, 4]);
    assert(dcomp(1, 3) > 0);
    expect_equal(dcomp(1, 2), dcomp(1, 1) / dcomp(1, 3));
    expect_equal(dcomp(1, 4), 0);

    expect_equal(numel(arr.regions), 1);
    region = arr.regions(1);
    expect_equal(region.name, 'solid');
    expect_equal(region.num_cells, 1);
    assert(region.components(1, 1) < dcomp(1, 1)); % one cell < both

    % [] names means every cell_data array (there's exactly one: material).
    expect_equal(numel(mio.data_integrate(m)), 1);

    % A point_data-only name fails, naming the fix.
    expect_error(@() mio.data_integrate(m, 'temperature'), 'point_data_to_cell_data');
end

% --- the settings pipeline -------------------------------------------------

function t_pipeline()
    % Behaviour follows the build: with the JSON parser a bad document fails
    % naming the offending op; without it every entry point fails naming
    % -DMESHIOPLUSPLUS_WITH_JSON=ON. Never a missing symbol either way.
    bad = ['{"Input": {"Path": "a"}, "Output": {"Path": "b"},', ...
           ' "Operations": [{"Op": "Nope"}]}'];
    if mio.pipeline_has_json()
        expect_error(@() mio.pipeline_run_json(bad), 'Nope');
        dir = tempname();
        mkdir(dir);
        cleanup = onCleanup(@() rmdir(dir, 's'));
        inp = fullfile(dir, 'in.vtu');
        out = fullfile(dir, 'out.vtu');
        mio.write(fixture(), inp);
        settings = fullfile(dir, 'settings.json');
        fid = fopen(settings, 'w');
        fprintf(fid, ['{"Input": {"Path": "%s"}, "Operations": [{"Op": "Quality"}], ' ...
                      '"Output": {"Path": "%s"}}\n'], inp, out);
        fclose(fid);
        mio.pipeline_run_file(settings);
        report = mio.pipeline_run_file_report(settings);
        assert(ischar(report) && size(report, 1) == 1);
        assert(~isempty(strfind(report, '"op":"Quality"')));
        assert(~isempty(strfind(report, '"warnings":[]')));
        expect_error(@() mio.pipeline_run_json_report(bad), 'Nope');
        back = mio.read(out);
        assert(any(strcmp('quality:scaled_jacobian', mio.cell_data_names(back))));
        delete(out);
        text = fileread(settings);
        mio.pipeline_run_json(text);
        assert(exist(out, 'file') == 2);
        report2 = mio.pipeline_run_json_report(text);
        assert(~isempty(strfind(report2, '"op":"Quality"')));
    else
        expect_error(@() mio.pipeline_run_json(bad), 'MESHIOPLUSPLUS_WITH_JSON');
        expect_error(@() mio.pipeline_run_json_report(bad), 'MESHIOPLUSPLUS_WITH_JSON');
    end
end

function t_pipeline_v2()
    skip_if_not(mio.pipeline_has_json(), 'library built without the JSON parser');
    dir = tempname();
    mkdir(dir);
    cleanup = onCleanup(@() rmdir(dir, 's'));
    input = fullfile(dir, 'in.vtu');
    mio.write(fixture(), input, 'codec', 'none');
    text = sprintf(['{"Version":2,"Input":{"Path":"%s"},"Operations":[{"Op":"Merge",' ...
                    '"Inputs":["%s"]},{"Op":"Partition","Nparts":2,"Method":"sfc"}],' ...
                    '"Output":{"Pattern":"%s/part_{part}.vtu","Codec":"none"}}'], ...
                   input, input, dir);
    report = mio.pipeline_run_json_report(text);
    assert(~isempty(strfind(report, '"op":"Merge"')));
    assert(~isempty(strfind(report, '"NumPieces":2')));
    assert(exist(fullfile(dir, 'part_0.vtu'), 'file') == 2);
    assert(exist(fullfile(dir, 'part_1.vtu'), 'file') == 2);
end

function t_sequence_report()
    bad = '{"Input": {"Path": "a"}, "Output": {"Path": "b"}, "Operations": [{"Op": "Nope"}]}';
    if ~mio.pipeline_has_json()
        expect_error(@() mio.sequence_pipeline_run_json_report(bad), 'MESHIOPLUSPLUS_WITH_JSON');
        return;
    end
    expect_error(@() mio.sequence_pipeline_run_json_report(bad), 'Nope');
    dir = tempname();
    mkdir(dir);
    cleanup = onCleanup(@() rmdir(dir, 's'));
    inp = fullfile(dir, 'in.vtu');
    out = fullfile(dir, 'out.vtu');
    mio.write(fixture(), inp);
    % A document using none of the sequence keys behaves as a plain pipeline.
    text = sprintf(['{"Input": {"Path": "%s"}, "Operations": [{"Op": "Quality"}], ' ...
                    '"Output": {"Path": "%s"}}'], inp, out);
    report = mio.sequence_pipeline_run_json_report(text);
    assert(ischar(report) && ~isempty(strfind(report, '"op":"Quality"')));
    assert(exist(out, 'file') == 2);
    settings = fullfile(dir, 'settings.json');
    fid = fopen(settings, 'w');
    fprintf(fid, '%s\n', text);
    fclose(fid);
    report2 = mio.sequence_pipeline_run_file_report(settings);
    assert(~isempty(strfind(report2, '"op":"Quality"')));
end

% --- grids and generators --------------------------------------------------

function t_grid()
    g = mio.grid([2 2 2]);
    expect_equal(mio.num_points(g), 27);
    expect_equal(mio.num_cell_blocks(g), 1);
    p = mio.points(g);
    expect_equal(max(p, [], 1), [2 2 2]);

    s = mio.grid([1 1 1], 'origin', [1 2 3], 'spacing', [0.5 0.25 2]);
    expect_equal(max(mio.points(s), [], 1), [1.5 2.25 5]);

    % An empty lattice is a legal request, not an error.
    e = mio.grid([0 0 0]);
    expect_equal(mio.num_points(e), 0);

    expect_error(@() mio.grid([100 100 100], 'max_cells', 1000));
    expect_error(@() mio.grid([2 2]), 'dims');
end

function t_generators()
    b = mio.box('dims', [2 2 2], 'cell_type', 'tetra', 'order', 2);
    expect_equal(mio.num_cell_blocks(b), 1);

    s = mio.uv_sphere('nu', 8, 'nv', 4, 'cell_type', 'quad');
    expect_equal(mio.num_cell_blocks(s), 2);

    r = mio.rectangle('dims', [3 2], 'cell_type', 'triangle', 'split', 'center');
    expect_equal(mio.num_points(r), 18);

    ico = mio.icosa_sphere('subdivisions', 1);
    expect_equal(mio.num_points(ico), 42);

    d = mio.disk('n', 3, 'cell_type', 'triangle');
    expect_equal(mio.num_cell_blocks(d), 1);

    cyl = mio.cylinder('nxy', 2, 'nz', 2, 'cell_type', 'wedge');
    expect_equal(mio.num_cell_blocks(cyl), 1);

    ref = mio.reference_cells({'tetra', 'wedge'});
    expect_equal(mio.num_cell_blocks(ref), 2);

    % The planar domains: closed-form point counts, and an annulus sector
    % keeps its seam apart where a full ring welds it.
    expect_point_count(mio.reference_triangle('n', 3), 10);
    expect_point_count(mio.ngon('sides', 5, 'rings', 2), 16);
    expect_point_count(mio.annulus('ntheta', 12, 'nr', 3), 48);
    expect_point_count(mio.annulus('ntheta', 6, 'nr', 2, 'angle', pi / 2), 21);
    expect_point_count(mio.lshape('n', 3), 40);
    expect_point_count(mio.plate_with_hole('nt', 3, 'nr', 2), 36);
    expect_point_count(mio.honeycomb('dims', [3 2]), 22);
    expect_point_count(mio.line('n', 4, 'stop', [0 0 2]), 5);
    expect_point_count(mio.circle('n', 6, 'order', 2), 12);
    expect_point_count(mio.helix('n', 10), 11);
    expect_point_count(mio.sphere_points('n', 7), 7);
    expect_point_count(mio.box_points('n', 7, 'seed', 3), 7);
    expect_point_count(mio.tetra_sphere('subdivisions', 1), 10);
    expect_point_count(mio.octa_sphere('subdivisions', 1), 18);
    expect_point_count(mio.cube_sphere('n', 2), 26);
    expect_point_count(mio.geo_sphere('n', 2), 42);
    expect_point_count(mio.torus('nu', 6, 'nv', 4), 24);
    expect_point_count(mio.moebius('nu', 6, 'nv', 2), 18);
    expect_point_count(mio.tube('ntheta', 6, 'nz', 2), 18);
    expect_point_count(mio.cone('ntheta', 6, 'nz', 2), 20);
    expect_point_count(mio.capsule('ntheta', 6, 'nz', 1, 'ncap', 2), 26);
    expect_point_count(mio.box_surface('dims', [1 2 3]), 24);
    expect_point_count(mio.platonic('dodecahedron'), 20);
    expect_point_count(mio.ball('n', 1), 27);
    expect_point_count(mio.spherical_shell('n', 2, 'nr', 1), 52);
    expect_point_count(mio.pipe('ntheta', 6, 'nr', 1, 'nz', 2), 36);
    expect_point_count(mio.fichera('n', 1), 26);
    expect_point_count(mio.solid_torus('n', 2, 'nu', 5), 45);
    expect_point_count(mio.box('order', 2, 'complete', true), 27);
    expect_point_count(mio.ball('n', 1, 'order', 3), 343);
    expect_error(@() mio.platonic('sphere'), 'platonic');
    expect_error(@() mio.plate_with_hole('radius', 0.9));

    expect_error(@() mio.box('cell_type', 'blob'), 'cell_type');
    expect_error(@() mio.rectangle('split', 'diagonal'), 'split');
    expect_error(@() mio.box('dims', [100 100 100], 'max_cells', 1000));

    % Positional options follow R's declaration order.
    expect_point_count(mio.rectangle([2 2], [2 1]), 9);
    pr = mio.points(mio.rectangle([2 2], [2 1], [1 1 0]));
    expect_equal(max(pr, [], 1), [3 2 0]);
end

function expect_point_count(m, n)
    expect_equal(mio.num_points(m), n);
end

function t_generator_smoke()
    names = {'rectangle', 'box', 'disk', 'uv_sphere', 'icosa_sphere', 'cylinder', ...
             'reference_triangle', 'ngon', 'annulus', 'lshape', 'plate_with_hole', ...
             'honeycomb', 'line', 'circle', 'helix', 'sphere_points', 'box_points', ...
             'tetra_sphere', 'octa_sphere', 'cube_sphere', 'geo_sphere', 'torus', ...
             'moebius', 'tube', 'cone', 'capsule', 'box_surface', 'platonic', 'ball', ...
             'spherical_shell', 'pipe', 'fichera', 'solid_torus'};
    for k = 1:numel(names)
        m = feval(['mio.' names{k}]);
        assert(isa(m, 'mio.MeshHandle'), names{k});
        assert(mio.num_points(m) > 0, names{k});
        assert(mio.num_cell_blocks(m) >= 1, names{k});
        p = mio.points(m);
        expect_equal(size(p, 1), mio.num_points(m));
    end
end

function t_reference_cells()
    a = mio.reference_cells('tetra,wedge');
    expect_equal(mio.cell_block_types(a), {'tetra'; 'wedge'});
    pairs = mio.reference_cells('types', 'tetra', 'pairs', true);
    expect_equal(mio.num_cells(pairs), 2);
    every = mio.reference_cells();
    assert(mio.num_cell_blocks(every) > 2);
end

% --- voxelize, octree, signed distance -------------------------------------

function t_voxelize()
    cube = cube_surface();

    q = mio.surface_watertight_check(cube);
    assert(q.watertight);
    expect_equal(q.boundary_edges, 0);
    expect_equal(q.non_manifold_edges, 0);

    v = mio.voxelize(cube, 'resolution', [4 4 4]);
    expect_equal(v.num_occupied, 64);
    expect_equal(v.dims, [4 4 4]);
    expect_equal(v.spacing(1), 0.25);
    expect_equal(size(v.origin), [1, 3]);

    inside = mio.voxelize(cube, 'resolution', [5 5 5], ...
        'bounds', [-0.5 -0.5 -0.5 1.5 1.5 1.5], 'fill', 'inside', 'watertight_check', 'off');
    expect_equal(inside.num_occupied, 27);

    occ = mio.voxelize(cube, 'cell_size', 0.5, 'attach_occupancy', true);
    assert(any(strcmp('voxel:occupancy', mio.cell_data_names(occ.mesh))));

    tree = mio.octree(cube, 'max_level', 3, 'bounds', [-0.5 -0.5 -0.5 1.5 1.5 1.5], ...
        'fill', 'inside', 'watertight_check', 'off');
    expect_equal(tree.dimension, 3);
    assert(tree.num_leaves > 0);
    expect_equal(tree.finest_level, 3);
    expect_equal(size(tree.origin), [1, 3]);
    expect_error(@() mio.octree(cube, 'balance', 'diagonal'), 'balance');
    expect_error(@() mio.octree(cube, 'transition', 'wedges'), 'transition');
    templ = mio.octree(cube, 'max_level', 3, 'bounds', [-0.5 -0.5 -0.5 1.5 1.5 1.5], ...
        'transition', 'templates');
    expect_equal(templ.num_hanging, 0);
    rooted = mio.octree(cube, 'root_cells', [2 2 2], 'max_level', 1);
    assert(rooted.num_leaves >= 8);

    % An unknown fill fails by name rather than silently defaulting.
    expect_error(@() mio.voxelize(cube, 'resolution', [2 2 2], 'fill', 'solid'), 'fill');
    flat = mio.voxelize(cube, 'resolution', [4 4 0], 'plane', 'xy');
    expect_equal(flat.dims, [4 4 1]);
    expect_error(@() mio.voxelize(cube, 'resolution', [4 4 0], 'plane', 'uv'), 'plane');
    expect_error(@() mio.voxelize(cube, 'resolution', [4 4]), 'resolution');
end

function t_octree_mesh()
    cube = cube_surface();
    % octree_mesh: graded tetrahedra inside the cube, then its castellated leaves.
    tets = mio.octree_mesh(cube, 'max_size', 0.5, 'min_size', 0.125);
    expect_equal(tets.output, 'tet');
    assert(any(strcmp('tetra', mio.cell_block_types(tets.mesh))));
    assert(tets.num_tetrahedra > 0);
    expect_equal(tets.num_inverted, 0);
    assert(iscell(tets.warnings));
    expect_equal(numel(tets.warnings), tets.num_warnings);
    expect_equal(size(tets.origin), [1, 3]);
    expect_equal(size(tets.root_cells), [1, 3]);
    assert(islogical(tets.guarantee_applies));
    cells = mio.octree_mesh(cube, 'output', 'castellated', 'max_size', 0.5, 'min_size', 0.125);
    expect_equal(cells.output, 'castellated');
    assert(any(strcmp('hexahedron', mio.cell_block_types(cells.mesh))));
    assert(cells.num_hexahedra > 0);
    % A max_level too coarse for min_size is a warning, one entry per line.
    coarse = mio.octree_mesh(cube, 'max_size', 0.5, 'min_size', 0.125, 'max_level', 1);
    expect_equal(coarse.num_warnings, 1);
    assert(~isempty(regexp(coarse.warnings{1}, 'coarsens min_size', 'once')));
    expect_error(@() mio.octree_mesh(cube, 'output', 'prisms'), 'output');
    expect_error(@() mio.octree_mesh(cube, 'fill', 'surface'), 'fill');
end

function t_distance()
    cube = cube_surface();

    % The centre is 0.5 inside; the other two are 1.0 outside.
    pts = [0.5 0.5 0.5; 2 0.5 0.5; -1 0.5 0.5];
    d = mio.sample_distance(cube, pts, 'watertight_check', 'off');
    expect_equal(size(d), [3, 1]);
    expect_equal(d(1), -0.5, 1e-12);
    expect_equal(d(2), 1.0, 1e-12);
    expect_equal(d(3), 1.0, 1e-12);
    u = mio.sample_distance(cube, pts, 'sign', 'unsigned', 'watertight_check', 'off');
    expect_equal(u(1), 0.5, 1e-12);
    w = mio.sample_distance(cube, pts, 'sign', 'winding-number', 'watertight_check', 'off');
    expect_equal(w(1), -0.5, 1e-12);
    b = mio.sample_distance(cube, pts, 'band', 0.25, 'watertight_check', 'off');
    assert(abs(b(2)) <= 0.25 + 1e-12);

    g = mio.grid([2 2 2], 'origin', [-0.5 -0.5 -0.5]);
    f = mio.distance_to_surface(g, cube, 'watertight_check', 'off', 'record_inside', true);
    assert(f.quality.watertight);
    expect_equal(f.num_banded, 0);
    assert(any(strcmp('sdf:distance', mio.point_data_names(f.mesh))));
    assert(any(strcmp('sdf:inside', mio.point_data_names(f.mesh))));

    expect_error(@() mio.sample_distance(cube, pts, 'sign', 'magic'), 'sign');
    expect_error(@() mio.sample_distance(cube, pts(:, 1:2)), 'n x 3');
end

function t_distance_center()
    cube = cube_surface();
    dom = mio.grid([4 4 4], 'origin', [-0.5 -0.5 -0.5], 'spacing', [0.5 0.5 0.5]);
    f = mio.distance_to_surface(dom, cube, 'location', 'center', 'watertight_check', 'off');
    dist = mio.cell_data(f.mesh, 'sdf:distance');
    expect_equal(size(dist), [64, 1]);
    expect_equal(sum(dist < 0), 8);
    expect_error(@() mio.distance_to_surface(dom, cube, 'location', 'edge'), 'location');
end

function t_compute_sdf()
    cube = cube_surface();

    sdf = mio.compute_sdf(cube, 'resolution', [4 4 4], 'watertight_check', 'off');
    expect_equal(sdf.dims, [4 4 4]);
    expect_equal(sdf.max_depth, 0);
    assert(any(strcmp('sdf:distance', mio.point_data_names(sdf.mesh))));
    info = mio.cell_block_info(sdf.mesh, 1);
    expect_equal(info.num_cells, 64);
    assert(isstruct(sdf.quality) && islogical(sdf.quality.watertight));

    % The octree refines only near the surface: more than the root, far less
    % than the uniform grid of the same finest resolution.
    tree = mio.compute_sdf(cube, 'structure', 'octree', 'root_resolution', 4, ...
        'max_depth', 2, 'watertight_check', 'off');
    expect_equal(tree.max_depth, 2);
    info = mio.cell_block_info(tree.mesh, 1);
    assert(info.num_cells > 64 && info.num_cells < 4096);

    % resolution/cell_size size a voxel grid; an octree's finest cell is
    % already determined by root_resolution and max_depth.
    expect_error(@() mio.compute_sdf(cube, 'structure', 'octree', 'resolution', [4 4 4]));
    expect_error(@() mio.compute_sdf(cube, 'structure', 'quadtree'), 'structure');
end

function t_remesh_volume()
    cube = cube_surface();
    % Unlike remesh, remesh_volume accepts the closed surface directly and
    % returns a tetra mesh.
    r = mio.remesh_volume(cube, 'cell_size', 0.4, 'watertight_check', 'off');
    expect_equal(mio.cell_block_type(r.mesh, 1), 'tetra');
    assert(r.num_tets > 0);
    expect_equal(r.num_tets, mio.num_cells(r.mesh));
    assert(r.num_vertices_warped >= 0);
    assert(r.num_tets_rejected >= 0);
    assert(r.num_non_manifold_edges >= 0);
    assert(r.input_quality.watertight);

    res = mio.remesh_volume(cube, 'resolution', [3 3 3], 'warp_fraction', 0, ...
        'watertight_check', 'off');
    assert(res.num_tets > 0);
    expect_equal(res.num_vertices_warped, 0);

    % Neither resolution nor cell_size given is rejected by name.
    expect_error(@() mio.remesh_volume(cube));
end

% --- select ----------------------------------------------------------------

function t_select()
    m = cube_surface();
    sel = mio.select(m, 'top', 'seed', [0.5 0.5 2]);
    expect_equal(sel.num_selected, 2);
    expect_equal(sel.num_grown, 2);
    rg = mio.regions(sel.mesh);
    expect_equal(rg(end).name, 'top');
    low = mio.select(m, 'low', 'kind', 'point', 'box', [-1 -1 -1 2 2 0.5]);
    expect_equal(low.num_selected, 4);
    expect_equal(low.num_grown, -1);
    one = mio.select(m, 'all', 'component', 1, 'quality_metric', 'min_angle', ...
        'quality_op', '>', 'quality_value', 0);
    expect_equal(one.num_selected, 12);
    expect_equal(one.num_components, 1);
    ball = mio.select(m, 'ball', 'kind', 'point', 'sphere', [0 0 0 0.1]);
    expect_equal(ball.num_selected, 1);
    half = mio.select(m, 'half', 'kind', 'point', 'plane', [0 0 0.5 0 0 1]);
    expect_equal(half.num_selected, 4);
    mio.add_point_data(m, 'z', mio.points(m) * [0; 0; 1]);
    hi = mio.select(m, 'hi', 'kind', 'point', 'where_array', 'z', 'where_op', '~=', ...
        'where_value', 0);
    expect_equal(hi.num_selected, 4);
    expect_error(@() mio.select(m, 'none'));
    expect_error(@() mio.select(m, 'x', 'kind', 'edge', 'component', 1), 'kind');
    expect_error(@() mio.select(m, 'x', 'box', [0 0 0]), 'box');
    expect_error(@() mio.select(m, 'x', 'component', 1, 'where_op', '<>'), 'where_op');
end

% --- provenance ------------------------------------------------------------

function t_provenance()
    % Outside a scope the recorders are no-ops.
    mio.provenance_note('dtype', 'nothing to record');
    mio.provenance_source('a.msh', 'gmsh');
    mio.provenance_target('vtu');

    m = fixture();
    path = [tempname() '.vtu'];
    cleanup = onCleanup(@() delete_if(path));
    v = mio.with_provenance(@() write_with_notes(m, path));
    expect_equal(v, 42);
    assert(exist(path, 'file') == 2);
    % The richer block carries the recorded source; a plain write does not.
    assert(~isempty(strfind(fileread(path), 'bracket.msh')));
    expect_equal(mio.points(mio.read(path)), fixture_points());
    plain = [tempname() '.vtu'];
    cleanup2 = onCleanup(@() delete_if(plain));
    mio.write(m, plain, 'codec', 'none');
    assert(isempty(strfind(fileread(plain), 'bracket.msh')));

    [a, b] = mio.with_provenance(@() deal(1, 2), 'mode', 'off');
    expect_equal([a, b], [1, 2]);
    % No output requested: the function still runs.
    mio.with_provenance(@() mio.provenance_note('x', 'y'));

    % The scope closes even when the function throws: a new scope opens after.
    expect_error(@() mio.with_provenance(@() error('custom:boom', 'boom')), 'boom');
    expect_equal(mio.with_provenance(@() 7), 7);
    expect_error(@() mio.with_provenance(@() 1, 'mode', 'loud'), 'mode');
    expect_error(@() mio.with_provenance(42), 'function handle');
end

function v = write_with_notes(m, path)
    mio.provenance_source('bracket.msh', 'gmsh');
    mio.provenance_note('regions-dropped', 'Side regions have no VTU equivalent');
    mio.write(m, path, 'codec', 'none');
    mio.provenance_target('vtu', 'encoding', 'binary', 'codec', 'none');
    v = 42;
end

% --- errors ----------------------------------------------------------------

function t_errors()
    m = fixture();
    expect_error(@() mio.data_rename(m, 'point', 'missing', 'x'));
    try
        mio.data_rename(m, 'point', 'missing', 'x');
    catch err
        % Octave prefixes the MEX name; the C library's message follows.
        assert(~isempty(strfind(err.message, mio.last_error())));
        assert(strncmp(err.identifier, 'meshioplusplus:', 15));
    end
end

function delete_if(path)
    if exist(path, 'file') == 2
        delete(path);
    end
end
