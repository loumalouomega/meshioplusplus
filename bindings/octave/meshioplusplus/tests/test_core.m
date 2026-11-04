function tests = test_core()
%TEST_CORE  Lifecycle, I/O, accessors and handle safety (mesh.c).
    tests = {
        'build information is available', @t_build_info
        'a mesh can be built and inspected', @t_inspect
        'rows are entities and values survive unchanged', @t_layout
        'data arrays keep their class', @t_types
        'a file round trip preserves everything', @t_roundtrip
        'struct conversion round-trips', @t_struct
        'released, stale and foreign handles are rejected', @t_handles
        'optional arguments are positional or name/value', @t_opts
        'ragged polygon and polyhedron blocks round-trip', @t_ragged
        'named regions round-trip with the facet rule', @t_regions
        'read_metadata summarises a file', @t_metadata
        'selective reads honour arrays and points_only', @t_selective
        'MED named meshes preserve file order', @t_med
        'Gmsh 4.0 reads natively with sparse node tags', @t_gmsh40
        'errors carry the C library message', @t_errors
    };
end

function t_build_info()
    assert(~isempty(mio.version()));
    assert(any(strcmp(mio.mesh_backend(), {'meshio', 'native', 'kratos'})));
    assert(mio.format_readable('vtu'));
    assert(mio.format_writable('vtu'));
    assert(mio.format_writable('openfoam'));
    assert(~mio.format_readable('nonexistent'));
    expect_equal(mio.cell_type_num_nodes('tetra10'), 10);
    expect_equal(mio.cell_type_dimension('triangle'), 2);
    expect_equal(mio.cell_type_num_nodes('polygon'), -1);
end

function t_inspect()
    m = fixture();
    expect_equal(mio.num_points(m), 5);
    expect_equal(mio.point_dim(m), 3);
    expect_equal(mio.num_cell_blocks(m), 1);
    expect_equal(mio.num_cells(m), 2);
    expect_equal(mio.cell_block_type(m, 1), 'tetra');
    expect_equal(mio.cell_block_types(m), {'tetra'});
    info = mio.cell_block_info(m, 1);
    expect_equal(info.num_cells, 2);
    expect_equal(info.nodes_per_cell, 4);
    assert(~info.is_ragged);
    expect_equal(mio.point_data_names(m), {'displacement'; 'temperature'});
    expect_equal(mio.cell_data_names(m), {'material'});
    expect_equal(mio.cell_data_num_blocks(m, 'material'), 1);
end

function t_layout()
    m = fixture();
    p = mio.points(m);
    expect_equal(size(p), [5, 3]);
    expect_equal(p, fixture_points());
    expect_equal(p(2, :), [1.1, 0.2, 0.3]);
    expect_equal(mio.connectivity(m), fixture_conn());
    expect_equal(mio.connectivity_raw(m), fixture_conn() - 1);
    expect_equal(mio.point_data(m, 'displacement'), fixture_vec());
    expect_equal(mio.point_data(m, 'temperature'), [1; 2; 3; 4; 5]);
    expect_equal(mio.cell_data(m, 'material'), [7; 9]);
end

function t_types()
    m = fixture();
    mio.add_point_data(m, 'ids', int32([5; 4; 3; 2; 1]));
    mio.add_point_data(m, 'flags', uint8([1; 0; 1; 0; 1]));
    mio.add_point_data(m, 'f32', single(fixture_vec()));
    mio.add_field_data(m, 'tensor', reshape(1:12, [2, 3, 2]));
    assert(isa(mio.point_data(m, 'ids'), 'int32'));
    assert(isa(mio.point_data(m, 'flags'), 'uint8'));
    assert(isa(mio.point_data(m, 'f32'), 'single'));
    expect_equal(mio.point_data(m, 'ids'), int32([5; 4; 3; 2; 1]));
    expect_equal(mio.point_data(m, 'f32'), single(fixture_vec()));
    expect_equal(mio.field_data(m, 'tensor'), reshape(1:12, [2, 3, 2]));
end

function t_roundtrip()
    m = fixture();
    path = [tempname() '.vtu'];
    cleanup = onCleanup(@() delete_if(path));
    mio.write(m, path, 'codec', 'none');
    r = mio.read(path);
    expect_equal(mio.points(r), fixture_points());
    expect_equal(mio.connectivity(r), fixture_conn());
    expect_equal(mio.point_data(r, 'displacement'), fixture_vec());
    expect_equal(mio.cell_data(r, 'material'), [7; 9]);
    out = [tempname() '.msh'];
    cleanup2 = onCleanup(@() delete_if(out));
    mio.convert(path, out);
    expect_equal(mio.points(mio.read(out)), fixture_points());
    expect_equal(mio.sniff_format(out), 'gmsh');
    % Every write encoding reaches the C API (ASCII once tripped over a macro
    % in the C header shadowing its enumerator).
    for enc = {'ascii', 'binary', 'raw_appended'}
        p2 = [tempname() '.vtu'];
        mio.write(m, p2, 'encoding', enc{1});
        expect_equal(mio.points(mio.read(p2)), fixture_points());
        delete(p2);
    end
end

function t_struct()
    m = fixture();
    s = m.toStruct();
    expect_equal(s.points, fixture_points());
    expect_equal(s.cells(1).type, 'tetra');
    expect_equal(s.cells(1).data, fixture_conn());
    expect_equal(s.point_data(1).name, 'displacement');
    m2 = mio.MeshHandle.fromStruct(s);
    s2 = m2.toStruct();
    assert(isequal(s, s2));
end

function t_handles()
    m = fixture();
    assert(mio.is_open(m));
    mio.release(m);
    assert(~mio.is_open(m));
    mio.release(m); % idempotent
    expect_error(@() mio.points(m), 'released');
    expect_error(@() mio.points(42), 'mio.MeshHandle');
    forged = mio.MeshHandle(uint64(123456789));
    expect_error(@() mio.points(forged), 'released');
    assert(~forged.is_open());
    % Freed by the destructor when the last reference goes.
    a = fixture();
    b = a;
    clear a;
    expect_equal(mio.num_points(b), 5);
end

function t_opts()
    m = fixture();
    expect_equal(mio.connectivity(m, 1), fixture_conn());
    expect_equal(mio.connectivity(m, 'block', 1), fixture_conn());
    expect_error(@() mio.connectivity(m, 'bogus', 1), 'unknown option');
    expect_error(@() mio.connectivity(m, 0), '1-based');
    expect_error(@() mio.write(m, 'x.vtu', 'encoding', 'nope'), 'encoding');
end

function t_ragged()
    m = mio.mesh();
    mio.set_points(m, [0 0 0; 1 0 0; 1 1 0; 0 1 0; 2 0 0; 2 1 0; 0.5 0.5 1]);
    mio.add_polygon_block(m, 'polygon', {[1 2 3 4], [2 5 6], [3 6 7 4 1]});
    rows = mio.polygon_block(m, 1);
    expect_equal(rows, {[1 2 3 4]; [2 5 6]; [3 6 7 4 1]});
    expect_error(@() mio.connectivity(m, 1), 'ragged');
    p = mio.mesh();
    mio.set_points(p, [0 0 0; 1 0 0; 0 1 0; 0 0 1]);
    mio.add_polyhedron_block(p, 'polyhedron', {{[1 3 2], [1 2 4], [2 3 4], [3 1 4]}});
    cells = mio.polyhedron_block(p, 1);
    expect_equal(numel(cells), 1);
    expect_equal(cells{1}{4}, [3 1 4]);
end

function t_regions()
    m = fixture();
    mio.add_region(m, 'inlet', 'point', [1 2 3], 'dim', 2, 'tag', 7);
    mio.add_region(m, 'walls', 'side', [1 0; 2 3]);
    r = mio.regions(m);
    names = {r.name};
    i = find(strcmp(names, 'inlet'));
    expect_equal(r(i).kind, 'point');
    expect_equal(r(i).entries, [1; 2; 3]);
    expect_equal(r(i).tag, 7);
    j = find(strcmp(names, 'walls'));
    expect_equal(r(j).entries, [1 0; 2 3]);
    u = mio.edit_regions(m, 'rename', {'inlet'}, 'output', 'inflow');
    ru = mio.regions(u);
    assert(any(strcmp({ru.name}, 'inflow')));
    mio.remove_region(m, i);
    expect_equal(numel(mio.regions(m)), numel(r) - 1);
end

function t_metadata()
    m = fixture();
    path = [tempname() '.vtu'];
    cleanup = onCleanup(@() delete_if(path));
    mio.write(m, path);
    s = mio.read_metadata(path);
    expect_equal(s.num_points, 5);
    expect_equal(s.num_cells, 2);
    expect_equal(s.cell_block_types, {'tetra'});
    expect_equal(s.point_data_names, {'displacement'; 'temperature'});
end

function t_selective()
    m = fixture();
    path = [tempname() '.vtu'];
    cleanup = onCleanup(@() delete_if(path));
    mio.write(m, path);
    r = mio.read(path, 'arrays', {'temperature'});
    expect_equal(mio.point_data_names(r), {'temperature'});
    r = mio.read(path, 'arrays', {});
    expect_equal(numel(mio.point_data_names(r)), 0);
    r = mio.read(path, 'points_only', true);
    expect_equal(mio.num_points(r), 5);
    expect_equal(numel(mio.point_data_names(r)), 0);
end

function t_med()
    skip_if_not(mio.format_writable('med'), 'MED needs HDF5');
    first = fixture();
    second = fixture();
    path = [tempname() '.med'];
    cleanup = onCleanup(@() delete_if(path));
    mio.med_write_multi(path, {first, second}, {'z_mesh', 'a_mesh'});
    expect_equal(mio.med_mesh_names(path), {'z_mesh'; 'a_mesh'});
    sel = mio.med_read_named(path, 'a_mesh');
    expect_equal(mio.points(sel), mio.points(second));
    expect_error(@() mio.med_read_named(path, 'missing'), 'no mesh named');
    expect_error(@() mio.med_write_multi(path, {first, second}, {'same', 'same'}));
end

function t_gmsh40()
    path = [tempname() '.msh'];
    cleanup = onCleanup(@() delete_if(path));
    fid = fopen(path, 'w');
    fprintf(fid, '%s\n', '$MeshFormat', '4.0 0 8', '$EndMeshFormat', '$Nodes', '1 3', ...
            '22 2 0 3', '30 1 0 0', '10 0 0 0', '20 0 1 0', '$EndNodes', '$Elements', ...
            '1 1', '22 2 2 1', '90 10 30 20', '$EndElements');
    fclose(fid);
    m = mio.read(path, 'gmsh');
    expect_equal(mio.num_points(m), 3);
    expect_equal(mio.cell_block_type(m, 1), 'triangle');
    expect_equal(mio.connectivity(m, 1), [2 1 3]);
end

function t_errors()
    expect_error(@() mio.read('/definitely/not/here.vtu'));
    assert(~isempty(mio.last_error()));
    try
        mio.read('/definitely/not/here.vtu');
    catch err
        assert(strncmp(err.identifier, 'meshioplusplus:', 15));
    end
end

function delete_if(path)
    if exist(path, 'file')
        delete(path);
    end
end
