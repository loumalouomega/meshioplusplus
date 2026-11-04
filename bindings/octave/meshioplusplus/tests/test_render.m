function tests = test_render()
%TEST_RENDER  Software rendering (render.c; v16.33.0, field rendering v16.34.0).
    tests = {
        'mio.render returns a frame with pixels, ids and notes', @t_frame
        'the frame is upright: row 1 is the top of the picture', @t_orientation
        'a mapped field reports its range, ticks and symmetric limits', @t_field
        'a cut-away clips the near face and caps the cut', @t_cutaway
        'streamlines follow a vector point array', @t_streamlines
        'text, PNG and file forms render', @t_forms
        'options are name/value pairs, a struct, or positional', @t_calling
        'bad options are refused by name', @t_bad_options
        'the terminal colour depth follows the environment variables', @t_color_depth
    };
end

function m = render_cube()
    % The R suite's cube: six quads and a point array 'x' (the x coordinate).
    pts = [0 0 0; 1 0 0; 1 1 0; 0 1 0; 0 0 1; 1 0 1; 1 1 1; 0 1 1];
    quads = [1 4 3 2; 5 6 7 8; 1 2 6 5; 4 8 7 3; 1 5 8 4; 2 3 7 6];
    m = mio.mesh();
    mio.set_points(m, pts);
    mio.add_cell_block(m, 'quad', quads);
    mio.add_point_data(m, 'x', pts(:, 1));
end

function t_frame()
    m = render_cube();
    f = mio.render(m, 'width', 48, 'height', 32, 'background', mio.rgba(255, 255, 255));
    expect_equal(f.width, 48);
    expect_equal(f.height, 32);
    assert(isa(f.rgba, 'uint8'));
    expect_equal(size(f.rgba), [32, 48, 4]);
    expect_equal(size(f.cell_ids), [32, 48]);
    expect_equal(double(squeeze(f.rgba(1, 1, :))), [255; 255; 255; 255]);
    assert(sum(f.cell_ids(:) > 0) > 100);
    % 1-based global cell ids, 0 where no cell is drawn.
    assert(all(f.cell_ids(:) >= 0 & f.cell_ids(:) <= mio.num_cells(m)));
    assert(f.cell_ids(1, 1) == 0);
    assert(isempty(f.range));
    assert(iscellstr(f.notes));
    a = mio.frame_array(f);
    expect_equal(size(a), [32, 48, 4]);
    expect_equal(a(1, 1, 1), 1);
    % Background pixels are exactly where no cell was drawn.
    bg = all(f.rgba == 255, 3);
    assert(isequal(bg(f.cell_ids > 0), false(nnz(f.cell_ids > 0), 1)));
end

function t_orientation()
    % A triangle wide at the top (y = 1) and pointed at the bottom (y = 0),
    % seen down +z: the upright picture covers more of its first rows than
    % of its last ones.
    m = mio.mesh();
    mio.set_points(m, [0 1 0; 1 1 0; 0.5 0 0]);
    mio.add_cell_block(m, 'triangle', [1 2 3]);
    f = mio.render(m, 'width', 40, 'height', 40, 'view', '+z', 'shading', 'none');
    covered = sum(f.cell_ids > 0, 2);
    rows = find(covered > 0);
    assert(covered(rows(1) + 1) > covered(rows(end) - 1));
    % The PNG holds the same pixels in the same orientation.
    png = mio.render_png(m, 'width', 40, 'height', 40, 'view', '+z', 'shading', 'none');
    path = [tempname() '.png'];
    cleanup = onCleanup(@() delete_if(path));
    fid = fopen(path, 'wb');
    fwrite(fid, png);
    fclose(fid);
    try
        [img, ~, alpha] = imread(path);
    catch
        return; % no image reader in this installation
    end
    expect_equal(img, f.rgba(:, :, 1:3));
    expect_equal(alpha, f.rgba(:, :, 4));
end

function t_field()
    m = render_cube();
    g = mio.render(m, 'width', 40, 'height', 40, 'color_by', 'x', 'colorbar', true, ...
                   'symmetric', true, 'clip', [5 95], 'view', '+z');
    assert(~isempty(g.range));
    expect_equal(size(g.range), [1, 2]);
    expect_equal(g.range(1), -g.range(2));
    assert(any(strncmp(g.notes, 'x: ', 3)));
    assert(any(strncmp(g.notes, 'ticks: ', 7)));
    % NaN leaves one end of the clip open.
    h = mio.render(m, 'width', 40, 'height', 40, 'color_by', 'x', 'clip', [NaN 95]);
    assert(~isempty(h.range));
end

function t_cutaway()
    m = render_cube();
    f = mio.render(m, 'width', 32, 'height', 32, 'view', '-x', 'shading', 'none', ...
                   'background', mio.rgba(0, 0, 0), 'cutaway', [0.5 0 0 1 0 0], ...
                   'cutaway_tint', mio.rgba(10, 200, 30));
    % The cube closes, so the cut is capped (no cell, the fill colour), not tinted.
    expect_equal(f.cell_ids(17, 17), 0);
    assert(f.rgba(17, 17, 2) ~= 200);
    two = mio.render(m, 'width', 32, 'height', 32, ...
                     'cutaway', {[0.5 0 0 1 0 0], [0 0.5 0 0 1 0]});
    expect_equal(size(two.rgba), [32, 32, 4]);
    flat = mio.render(m, 'width', 32, 'height', 32, ...
                      'cutaway', [0.5 0 0 1 0 0, 0 0.5 0 0 1 0]);
    expect_equal(flat.rgba, two.rgba);
    expect_error(@() mio.render(m, 'cutaway', [1 2 3]), 'cutaway');
    expect_error(@() mio.render(m, 'cutaway', zeros(1, 6)), 'normal');
end

function t_streamlines()
    m = render_cube();
    pts = mio.points(m);
    mio.add_point_data(m, 'flow', [-(pts(:, 2) - 0.5), pts(:, 1) - 0.5, zeros(8, 1)]);
    opts = {'width', 64, 'height', 64, 'view', '+z', 'shading', 'none', ...
            'background', mio.rgba(0, 0, 0)};
    expect_equal(red(mio.render(m, opts{:})), 0);
    lines = mio.render(m, opts{:}, 'streamlines', 'flow', 'stream_color', mio.rgba(255, 0, 0));
    assert(red(lines) > 20);
    expect_error(@() mio.render(m, 'streamlines', 'nope'), 'nope');
    expect_error(@() mio.render(m, 'streamlines', 'flow', 'stream_seeds', 0), 'seeds');
end

function n = red(f)
    n = sum(sum(f.rgba(:, :, 1) == 255 & f.rgba(:, :, 2) == 0 & f.rgba(:, :, 3) == 0));
end

function t_forms()
    m = render_cube();
    txt = mio.render_text(m, 'cols', 20, 'rows', 6, 'color_depth', 'mono');
    assert(ischar(txt) && size(txt, 1) == 1);
    expect_equal(sum(txt == sprintf('\n')), 6);
    assert(~any(txt == char(27)));
    html = mio.render_text(m, 'cols', 20, 'rows', 6, 'format', 'html');
    assert(strncmp(html, '<!DOCTYPE html>', 15));
    plain = mio.render_text(m, 'cols', 20, 'rows', 6, 'format', 'plain', 'notes', false, ...
                            'encoding', 'ascii');
    assert(~any(plain == char(27)));

    png = mio.render_png(m, 'width', 30, 'height', 20);
    assert(isa(png, 'uint8'));
    expect_equal(size(png, 2), 1);
    expect_equal(png(1:8), uint8([137; 80; 78; 71; 13; 10; 26; 10]));
    expect_error(@() mio.render_png(m, 'compress', 10));

    path = [tempname() '.png'];
    cleanup = onCleanup(@() delete_if(path));
    mio.write_snapshot(path, m, 'width', 32, 'height', 24);
    expect_equal(read_bytes(path, 4), uint8([137; 80; 78; 71]));
    txtpath = [tempname() '.txt'];
    cleanup2 = onCleanup(@() delete_if(txtpath));
    mio.write_snapshot(txtpath, m, 'cols', 20, 'rows', 6, 'color_depth', 'mono');
    assert(numel(read_bytes(txtpath, Inf)) > 20);
    castpath = [tempname() '.cast'];
    cleanup3 = onCleanup(@() delete_if(castpath));
    mio.write_snapshot(castpath, m, 'cols', 16, 'rows', 4, 'cast_frames', 3, 'cast_fps', 6, ...
                       'cast_degrees', 90, 'png_compress', 0);
    head = char(read_bytes(castpath, 1)');
    expect_equal(head, '{');
    expect_error(@() mio.write_snapshot([tempname() '.vtu'], m));
end

function t_calling()
    m = render_cube();
    % Positional text options, as R's mio_render_text(m, 20L, 6L, ...).
    a = mio.render_text(m, 20, 6, 'halfblock', 'mono');
    b = mio.render_text(m, 'cols', 20, 'rows', 6, 'encoding', 'halfblock', 'color_depth', 'mono');
    expect_equal(a, b);
    % Positional text options, then render options.
    c = mio.render_text(m, 20, 6, 'view', '+z');
    d = mio.render_text(m, 'cols', 20, 'rows', 6, 'view', '+z');
    expect_equal(c, d);
    expect_error(@() mio.render_text(m, 20, 6, 'bogus', 1));
    expect_error(@() mio.render_text(m, 'bogus'), 'bogus');
    % A struct of options is the same as the pairs.
    s.width = 24;
    s.height = 16;
    s.view = '+x';
    f1 = mio.render(m, s);
    f2 = mio.render(m, 'width', 24, 'height', 16, 'view', '+x');
    expect_equal(f1.rgba, f2.rgba);
    % [] keeps an option at its default.
    f3 = mio.render(m, 'width', 24, 'height', 16, 'view', '+x', 'color_by', []);
    expect_equal(f3.rgba, f2.rgba);
    expect_equal(mio.render_png(m, 0, 'width', 8, 'height', 8), ...
                 mio.render_png(m, 'compress', 0, 'width', 8, 'height', 8));
    expect_error(@() mio.render(m, 1), 'name/value');
    expect_error(@() mio.render(m, 'width'), 'name/value');
    expect_error(@() mio.render(42), 'mio.MeshHandle');
    expect_error(@() mio.frame_array(42), 'mio.render');
    % Colours as [r g b] / [r g b a] and as the packed number.
    g1 = mio.render(m, 'width', 8, 'height', 8, 'background', [10 20 30]);
    g2 = mio.render(m, 'width', 8, 'height', 8, 'background', mio.rgba(10, 20, 30, 255));
    expect_equal(g1.rgba, g2.rgba);
    expect_equal(double(squeeze(g1.rgba(1, 1, :))), [10; 20; 30; 255]);
    expect_error(@() mio.render(m, 'background', [1 2]), 'background');
    expect_error(@() mio.render(m, 'background', [300 0 0]), '0..255');
    expect_error(@() mio.render(m, 'background', -1), '2\^32');
    % The remaining structured options are accepted.
    mio.render(m, 'width', 16, 'height', 16, 'pan', [0.1 0], 'light_dir', [0 0 1], ...
               'projection', 'perspective', 'edges', 'feature', 'color_by', 'x', ...
               'scale', 'linear', 'vmin', 0, 'vmax', 1, 'isolines', 0, 'iso_levels', [0.25 0.5]);
    mio.render(m, 'width', 16, 'height', 16, 'diagnostic', 'inverted');
    expect_error(@() mio.render(m, 'pan', 1), 'pan');
    expect_error(@() mio.render(m, 'projection', 'fisheye'), 'projection');
end

function t_bad_options()
    m = render_cube();
    expect_error(@() mio.render(m, 'color_by', 'nope'), 'nope');
    expect_error(@() mio.render(m, 'bogus', 1), 'unknown option ''bogus''');
    expect_error(@() mio.render(m, 'shading', 'plastic'), 'shading');
    expect_error(@() mio.render(m, 'width', 0), 'width and height');
    expect_error(@() mio.render_text(m, 'encoding', 'morse'), 'encoding');
end

function t_color_depth()
    expect_equal(mio.detect_color_depth('no_color', '1'), 'mono');
    expect_equal(mio.detect_color_depth('color_term', 'truecolor'), 'truecolor');
    expect_equal(mio.detect_color_depth('term', 'xterm-256color'), '256');
    expect_equal(mio.detect_color_depth('term', 'xterm'), '16');
    expect_equal(mio.detect_color_depth([], [], 'xterm'), '16');
    d = mio.detect_color_depth();
    assert(any(strcmp(d, {'truecolor', '256', '16', 'mono'})));
    expect_equal(mio.rgba(200, 197, 189), hex2dec('C8C5BDFF'));
    expect_equal(mio.rgba(1, 2, 3, 4), ((1 * 256 + 2) * 256 + 3) * 256 + 4);
end

function b = read_bytes(path, n)
    fid = fopen(path, 'rb');
    b = uint8(fread(fid, n, 'uint8'));
    fclose(fid);
end

function delete_if(path)
    if exist(path, 'file')
        delete(path);
    end
end
