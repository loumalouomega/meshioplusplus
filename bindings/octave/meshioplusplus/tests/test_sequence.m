function tests = test_sequence()
%TEST_SEQUENCE  Multi-file / transient sequences (sequence.c). Ports the
%   sequence-resampling test of test-meshioplusplus.R and adds smoke tests
%   for every other sequence entry point.
    tests = {
        'a glob opens a natural-sorted plan with times', @t_open_pattern
        'an explicit path list keeps its order unless sorted', @t_open_list
        'sequence resampling blends steps', @t_resample
        'fan-in and fan-out through a time series', @t_timeseries
        'a sequence settings document runs', @t_pipeline
        'sequence handles reject misuse and free idempotently', @t_handles
    };
end

% Steps step_0.vtu, step_1.vtu, ... step_<n-1>.vtu with point data u = 10 k.
function dir = write_steps(names)
    dir = tempname();
    mkdir(dir);
    for k = 1:numel(names)
        mk = fixture();
        mio.add_point_data(mk, 'u', repmat((k - 1) * 10, mio.num_points(mk), 1));
        mio.write(mk, fullfile(dir, names{k}));
        mio.release(mk);
    end
end

function t_open_pattern()
    % step_10 sorts after step_9: natural-numeric, not lexicographic.
    dir = write_steps({'step_0.vtu', 'step_9.vtu', 'step_10.vtu'});
    cleanup = onCleanup(@() rmdir(dir, 's'));
    seq = mio.sequence(fullfile(dir, 'step_*.vtu'), 'time_from', 'index');
    assert(isa(seq, 'mio.SequenceHandle'));
    expect_equal(mio.sequence_count(seq), 3);
    [~, name] = fileparts(mio.sequence_path(seq, 2));
    expect_equal(name, 'step_9');
    [~, name] = fileparts(mio.sequence_path(seq, 3));
    expect_equal(name, 'step_10');
    expect_equal(mio.sequence_step(seq, 1), 0);
    expect_equal(mio.sequence_time(seq, 3), 2);
    expect_equal(mio.sequence_time_source(seq, 1), 'index');
    m = mio.sequence_read(seq, 3); % step_10.vtu, the third file written
    expect_equal(mio.point_data(m, 'u'), repmat(20, 5, 1));
    expect_equal(mio.points(m), fixture_points());
    mio.release(m);
    expect_error(@() mio.sequence_read(seq, 4));
    mio.sequence_free(seq);

    % Explicit times, and the filename as the time source.
    seq = mio.sequence(fullfile(dir, 'step_*.vtu'), [], [0.5, 1.5, 2.5]);
    expect_equal(mio.sequence_time(seq, 2), 1.5);
    expect_equal(mio.sequence_time_source(seq, 2), 'explicit');
    mio.sequence_free(seq);
    seq = mio.sequence(fullfile(dir, 'step_*.vtu'), 'format', 'vtu', 'time_from', 'filename');
    expect_equal(mio.sequence_time(seq, 3), 10);
    expect_equal(mio.sequence_time_source(seq, 3), 'filename');
    mio.sequence_free(seq);
    expect_error(@() mio.sequence(fullfile(dir, 'step_*.vtu'), 'time_from', 'bogus'), 'time_from');
    expect_error(@() mio.sequence(fullfile(dir, 'none_*.vtu')));
    expect_error(@() mio.sequence(fullfile(dir, 'step_*.vtu'), 'times', [1, 2]));
end

function t_open_list()
    dir = write_steps({'a_2.vtu', 'a_1.vtu'});
    cleanup = onCleanup(@() rmdir(dir, 's'));
    paths = {fullfile(dir, 'a_2.vtu'), fullfile(dir, 'a_1.vtu')};
    seq = mio.sequence_list(paths, 'time_from', 'index');
    expect_equal(mio.sequence_count(seq), 2);
    expect_equal(mio.sequence_path(seq, 1), paths{1});
    mio.sequence_free(seq);
    seq = mio.sequence_list(paths, [], [], 'index', true);
    expect_equal(mio.sequence_path(seq, 1), paths{2});
    mio.sequence_free(seq);
    seq = mio.sequence_list(paths{1}, 'times', 7);
    expect_equal(mio.sequence_count(seq), 1);
    expect_equal(mio.sequence_time(seq, 1), 7);
    mio.sequence_free(seq);
    expect_error(@() mio.sequence_list({}), 'non-empty');
end

function t_resample()
    dir = write_steps({'step_0.vtu', 'step_1.vtu'});
    cleanup = onCleanup(@() rmdir(dir, 's'));
    n = 5;
    seq = mio.sequence(fullfile(dir, 'step_*.vtu'), 'time_from', 'index');
    mio.sequence_resample(seq, fullfile(dir, 'out_{index}.vtu'), [0, 0.5, 1]);
    r = mio.read(fullfile(dir, 'out_1.vtu'));
    expect_equal(double(mio.point_data(r, 'u')), repmat(5, n, 1));
    mio.release(r);
    expect_error(@() mio.sequence_resample(seq, fullfile(dir, 'bad_{index}.vtu'), 2));
    expect_error(@() mio.sequence_resample(seq, fullfile(dir, 'x_{index}.vtu'), 0, ...
        'method', 'cubic'), 'unknown resample method');
    % clamp takes the end step instead of failing; nearest picks a step.
    mio.sequence_resample(seq, fullfile(dir, 'c_{index}.vtu'), 2, 'clamp', true);
    r = mio.read(fullfile(dir, 'c_0.vtu'));
    expect_equal(double(mio.point_data(r, 'u')), repmat(10, n, 1));
    mio.release(r);
    mio.sequence_resample(seq, fullfile(dir, 'n_{index}.vtu'), 0.4, 'nearest', false, true, 'vtu');
    r = mio.read(fullfile(dir, 'n_0.vtu'));
    expect_equal(double(mio.point_data(r, 'u')), zeros(n, 1));
    mio.release(r);
    mio.sequence_free(seq);
end

function t_timeseries()
    dir = write_steps({'step_0.vtu', 'step_1.vtu', 'step_2.vtu'});
    cleanup = onCleanup(@() rmdir(dir, 's'));
    seq = mio.sequence(fullfile(dir, 'step_*.vtu'), 'time_from', 'index');
    series = fullfile(dir, 'series.xdmf');
    % ascii: XDMF's "XML" data format, so no HDF5 is needed.
    mio.sequence_to_timeseries(seq, series, 'ascii', true);
    mio.sequence_free(seq);
    meta = mio.read_metadata(series);
    expect_equal(numel(meta.time_values), 3);
    mio.timeseries_to_sequence(series, fullfile(dir, 'back_{step}.vtu'));
    back = mio.sequence(fullfile(dir, 'back_*.vtu'), 'time_from', 'index');
    expect_equal(mio.sequence_count(back), 3);
    m = mio.sequence_read(back, 3);
    expect_equal(double(mio.point_data(m, 'u')), repmat(20, 5, 1));
    mio.release(m);
    mio.sequence_free(back);
    % The default HDF data format needs HDF5; a build without it fails, never
    % writes a truncated series.
    seq = mio.sequence(fullfile(dir, 'step_*.vtu'));
    try
        mio.sequence_to_timeseries(seq, fullfile(dir, 'hdf.xdmf'), 'xdmf');
    catch err
        assert(~isempty(err.message));
    end
    % A format that cannot hold a series fails naming {step}.
    expect_error(@() mio.sequence_to_timeseries(seq, fullfile(dir, 'one.stl')));
    mio.sequence_free(seq);
    expect_error(@() mio.timeseries_to_sequence(series, fullfile(dir, 'nopattern.vtu')));
end

function t_pipeline()
    dir = write_steps({'step_0.vtu', 'step_1.vtu'});
    cleanup = onCleanup(@() rmdir(dir, 's'));
    pattern = strrep(fullfile(dir, 'step_*.vtu'), '\', '/');
    target = strrep(fullfile(dir, 'out_{step}.vtu'), '\', '/');
    json = sprintf('{"Version": 1, "Input": {"Pattern": "%s"}, "Output": {"Path": "%s"}}', ...
        pattern, target);
    try
        mio.sequence_pipeline_run_json(json);
    catch err
        skip_if_not(isempty(strfind(err.message, 'MESHIOPLUSPLUS_WITH_JSON')), ...
                    'requires JSON pipeline support');
        rethrow(err);
    end
    % {step} is zero-padded in a settings document's output.
    assert(exist(fullfile(dir, 'out_0000.vtu'), 'file') == 2);
    assert(exist(fullfile(dir, 'out_0001.vtu'), 'file') == 2);
    settings = fullfile(dir, 'settings.json');
    fid = fopen(settings, 'w');
    fprintf(fid, '%s', strrep(json, 'out_{step}', 'file_{step}'));
    fclose(fid);
    mio.sequence_pipeline_run_file(settings);
    assert(exist(fullfile(dir, 'file_0001.vtu'), 'file') == 2);
    expect_error(@() mio.sequence_pipeline_run_json('{"Version": 1, "Bogus": 1}'));
    expect_error(@() mio.sequence_pipeline_run_file(fullfile(dir, 'missing.json')));
end

function t_handles()
    dir = write_steps({'step_0.vtu'});
    cleanup = onCleanup(@() rmdir(dir, 's'));
    seq = mio.sequence(fullfile(dir, 'step_*.vtu'));
    assert(seq.is_open());
    m = fixture();
    % A mesh is not a sequence, and a sequence is not a mesh.
    expect_error(@() mio.sequence_count(m), 'mio.SequenceHandle');
    expect_error(@() mio.num_points(seq), 'mio.MeshHandle');
    % Free is idempotent and, as in R, a no-op for anything else.
    mio.sequence_free(m);
    assert(m.is_open());
    mio.sequence_free(42);
    mio.sequence_free(seq);
    mio.sequence_free(seq);
    assert(~seq.is_open());
    expect_error(@() mio.sequence_count(seq), 'released');
    mio.release(m);
end
