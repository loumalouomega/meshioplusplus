function write_snapshot(path, mesh, varargin)
%WRITE_SNAPSHOT  Render a mesh to a file chosen by extension.
%   mio.write_snapshot(PATH, MESH, 'name', value, ...) writes .png, .txt,
%   .ansi, .html or .cast (an asciinema orbit). Text forms are sized 'cols'
%   x 'rows'. Every option is a name/value pair (or one struct):
%     'cols', 'rows', 'encoding', 'color_depth', 'format', 'cell_aspect',
%     'tmux', 'notes'   as for mio.render_text (80, 24, 'halfblock',
%                       'truecolor', 'ansi', 2, false, true)
%     'png_compress'    PNG compression level (0)
%     'cast_frames'     frame count of the asciicast orbit (36)
%     'cast_fps'        its frame rate (12)
%     'cast_degrees'    its sweep in degrees (360)
%   Any other name/value pair is a render option, see mio.render.

    [o, pairs] = mio.internal.render_split(varargin, {'cols', 80, 'rows', 24, ...
        'encoding', 'halfblock', 'color_depth', 'truecolor', 'format', 'ansi', ...
        'cell_aspect', 2, 'tmux', false, 'notes', true, 'png_compress', 0, ...
        'cast_frames', 36, 'cast_fps', 12, 'cast_degrees', 360}, 'mio.write_snapshot', false);
    meshioplusplus_mex('mio_write_snapshot', mio.internal.str(path, 'path'), mesh, pairs, ...
        mio.internal.render_text_opts(o), double(o.png_compress), double(o.cast_frames), ...
        double(o.cast_fps), double(o.cast_degrees));
end
