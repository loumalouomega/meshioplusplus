function tests = test_series()
%TEST_SERIES  Transient XDMF and Exodus series writers (xdmf_series.c,
%   exodus_series.c).
    tests = {
        'a transient XDMF series round-trips', @t_xdmf_roundtrip
        'XDMF series flush, auto_flush and append', @t_xdmf_flush_append
        'XDMF series and meshes reject each other', @t_xdmf_handles
        'an XDMF series is finalized when its last reference goes', @t_xdmf_destructor
        'Exodus series own their grid and expose an idempotent lifecycle', @t_exodus
        'Exodus series reject foreign handles', @t_exodus_handles
    };
end

function t_xdmf_roundtrip()
    % The one writer mio.write cannot express: the grid goes out once and each
    % step is appended. 'XML' keeps everything in the single .xdmf, so this runs
    % against a library built without HDF5 too.
    dir = tempname();
    mkdir(dir);
    cleanup = onCleanup(@() rmdir_if(dir));
    path = fullfile(dir, 'series.xdmf');

    m = fixture();
    mio.add_region(m, 'anchors', 'point', [1; 5], 'dim', 0, 'tag', 7);
    mio.add_region(m, 'empty', 'cell', zeros(0, 1), 'dim', 3, 'tag', 21);
    mio.add_region(m, 'wall', 'side', [1, 1], 'dim', 2, 'tag', 9);
    expected_regions = mio.regions(m);
    s = mio.xdmf_series(path, 'data_format', 'XML');
    assert(isa(s, 'mio.XdmfSeries'));
    assert(mio.xdmf_series_is_open(s));
    expect_equal(mio.xdmf_series_num_steps(s), 0);
    assert(~mio.xdmf_series_finalized(s));

    mio.xdmf_series_write_points_cells(s, m);

    % temperature = t + node id, so no two steps share a value and a step
    % mix-up cannot pass.
    times = [0, 0.5, 1];
    step = mio.mesh();
    mio.set_points(step, fixture_points());
    mio.add_cell_block(step, 'tetra', fixture_conn());
    for t = times
        mio.add_point_data(step, 'temperature', t + (1:5)');
        mio.xdmf_series_write_data(s, t, step);
    end
    expect_equal(mio.xdmf_series_num_steps(s), 3);
    assert(~isempty(strfind(evalc('disp(s)'), '3 step')));

    % The .xdmf is buffered until finalize; nothing is readable before it.
    mio.xdmf_series_finalize(s);
    mio.xdmf_series_finalize(s); % idempotent
    assert(mio.xdmf_series_finalized(s));
    assert(~isempty(strfind(evalc('disp(s)'), 'finalized')));
    mio.xdmf_series_release(s);
    assert(~mio.xdmf_series_is_open(s));
    mio.xdmf_series_release(s); % idempotent
    expect_error(@() mio.xdmf_series_num_steps(s), 'released');
    assert(~isempty(strfind(evalc('disp(s)'), 'released')));

    meta = mio.read_metadata(path);
    expect_equal(meta.time_values, times');

    for k = 1:numel(times)
        back = mio.read(path, 'time_step', k - 1);
        expect_equal(mio.num_points(back), 5);
        expect_equal(mio.point_data(back, 'temperature'), times(k) + (1:5)');
        expect_equal(mio.regions(back), expected_regions);
    end

    % An unknown data format is an error carrying the C API's own message.
    expect_error(@() mio.xdmf_series(fullfile(dir, 'bad.xdmf'), 'data_format', 'NoSuchFormat'));
    expect_error(@() mio.xdmf_series(fullfile(dir, 'bad.xdmf'), 'XML', -1, 'sideways'), 'mode');
end

function t_xdmf_flush_append()
    dir = tempname();
    mkdir(dir);
    cleanup = onCleanup(@() rmdir_if(dir));
    path = fullfile(dir, 'flush.xdmf');
    m = fixture();

    % Positional options: data_format, gzip_level, mode, auto_flush.
    s = mio.xdmf_series(path, 'XML', -1, 'truncate', false);
    mio.xdmf_series_write_points_cells(s, m);
    mio.xdmf_series_write_data(s, 0, m);
    mio.xdmf_series_flush(s);
    % A flush publishes the steps so far while the series stays open.
    assert(~mio.xdmf_series_finalized(s));
    expect_equal(time_values(path), 0);
    mio.xdmf_series_write_data(s, 1, m);
    mio.xdmf_series_release(s);
    expect_equal(time_values(path), [0; 1]);

    % 'append' continues the series already at the path.
    s = mio.xdmf_series(path, 'data_format', 'XML', 'mode', 'append', 'auto_flush', true);
    expect_equal(mio.xdmf_series_num_steps(s), 2);
    mio.xdmf_series_write_data(s, 2, m);
    % auto_flush publishes the step without an explicit flush.
    expect_equal(time_values(path), [0; 1; 2]);
    s.release();
    assert(~mio.xdmf_series_is_open(s));
end

function t_xdmf_handles()
    % A mesh and a series are separately registered: neither is accepted for
    % the other, and a foreign value is rejected rather than dereferenced.
    dir = tempname();
    mkdir(dir);
    cleanup = onCleanup(@() rmdir_if(dir));
    m = fixture();
    expect_error(@() mio.xdmf_series_num_steps(m), 'mio.XdmfSeries');
    other = mio.xdmf_series(fullfile(dir, 'b.xdmf'), 'data_format', 'XML');
    expect_error(@() mio.num_points(other), 'mio.MeshHandle');
    expect_error(@() mio.xdmf_series_num_steps({}), 'mio.XdmfSeries');
    expect_error(@() mio.xdmf_series_release(m), 'mio.XdmfSeries');
    assert(~mio.xdmf_series_is_open(m));
    assert(~mio.xdmf_series_is_open(42));
    forged = mio.XdmfSeries(uint64(123456789));
    expect_error(@() mio.xdmf_series_num_steps(forged), 'released');
    assert(~mio.xdmf_series_is_open(forged));
    % Released before the tempdir goes away: a destructor running afterwards
    % would try to write the .xdmf into a directory that no longer exists.
    mio.xdmf_series_release(other);
end

function t_xdmf_destructor()
    dir = tempname();
    mkdir(dir);
    cleanup = onCleanup(@() rmdir_if(dir));
    path = fullfile(dir, 'gc.xdmf');
    m = fixture();
    s = mio.xdmf_series(path, 'data_format', 'XML');
    mio.xdmf_series_write_points_cells(s, m);
    mio.xdmf_series_write_data(s, 0.25, m);
    clear s;
    expect_equal(time_values(path), 0.25);
end

function t_exodus()
    path = [tempname() '.e'];
    cleanup = onCleanup(@() delete_if(path));
    if ~mio.format_writable('exodus')
        expect_error(@() mio.exodus_series(path), 'MESHIOPLUSPLUS_WITH_NETCDF');
        return;
    end
    m = fixture();
    s = mio.exodus_series(path);
    assert(isa(s, 'mio.ExodusSeries'));
    assert(mio.exodus_series_is_open(s));
    expect_equal(mio.exodus_series_num_steps(s), 0);
    assert(~mio.exodus_series_finalized(s));
    expect_error(@() mio.exodus_series_write_data(s, 0, m), 'write_points_cells');
    mio.exodus_series_write_points_cells(s, m);
    mio.exodus_series_write_data(s, 0.123456789012345, m);
    mio.exodus_series_write_data(s, 1.5, m);
    expect_equal(mio.exodus_series_num_steps(s), 2);
    assert(~isempty(strfind(evalc('disp(s)'), '2 step')));
    mio.exodus_series_flush(s);
    mio.exodus_series_finalize(s);
    mio.exodus_series_finalize(s);
    assert(mio.exodus_series_finalized(s));
    expect_error(@() mio.exodus_series_write_data(s, 2, m), 'open series');
    back = mio.read(path, 'time_step', -1);
    expect_equal(mio.point_data(back, 'temperature'), [1; 2; 3; 4; 5]);
    mio.exodus_series_release(s);
    mio.exodus_series_release(s);
    assert(~mio.exodus_series_is_open(s));
    expect_error(@() mio.exodus_series_num_steps(s), 'released');
    assert(~isempty(strfind(evalc('disp(s)'), 'released')));
end

function t_exodus_handles()
    m = fixture();
    expect_error(@() mio.exodus_series_num_steps(m), 'mio.ExodusSeries');
    expect_error(@() mio.exodus_series_release(m), 'mio.ExodusSeries');
    expect_error(@() mio.exodus_series_flush([]), 'mio.ExodusSeries');
    assert(~mio.exodus_series_is_open(m));
    assert(~mio.exodus_series_is_open('x'));
    forged = mio.ExodusSeries(uint64(987654321));
    expect_error(@() mio.exodus_series_finalize(forged), 'released');
    % An XDMF series is not an Exodus series, nor the other way round.
    dir = tempname();
    mkdir(dir);
    cleanup = onCleanup(@() rmdir_if(dir));
    x = mio.xdmf_series(fullfile(dir, 'x.xdmf'), 'data_format', 'XML');
    expect_error(@() mio.exodus_series_num_steps(x), 'mio.ExodusSeries');
    assert(~mio.exodus_series_is_open(x));
    mio.xdmf_series_release(x);
    if mio.format_writable('exodus')
        path = [tempname() '.e'];
        cleanup2 = onCleanup(@() delete_if(path));
        e = mio.exodus_series(path);
        expect_error(@() mio.xdmf_series_num_steps(e), 'mio.XdmfSeries');
        mio.exodus_series_release(e);
    end
end

function delete_if(path)
    if exist(path, 'file')
        delete(path);
    end
end

function rmdir_if(dir)
    if exist(dir, 'dir')
        delete(fullfile(dir, '*'));
        rmdir(dir);
    end
end

function t = time_values(path)
    meta = mio.read_metadata(path);
    t = meta.time_values;
end
