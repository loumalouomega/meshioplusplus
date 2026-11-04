function txt = render_text(mesh, varargin)
%RENDER_TEXT  Render a mesh as text for a terminal.
%   TXT = mio.render_text(MESH, ...) returns a char row (UTF-8 block
%   characters, or a graphics-protocol escape), ready for fprintf('%s', TXT).
%   Options, positional in this order or as name/value pairs:
%     'cols', 'rows'   size in terminal cells (80, 24)
%     'encoding'       'halfblock' (default), 'quadrant', 'sextant',
%                      'braille', 'ascii', or a graphics protocol: 'kitty',
%                      'iterm2' or 'sixel'
%     'color_depth'    'truecolor' (default), '256', '16' or 'mono'
%     'format'         'ansi' (SGR colour escapes; default), 'plain' (glyphs
%                      only) or 'html'
%     'cell_aspect'    height over width of a terminal cell (2)
%     'tmux'           wrap a graphics protocol for tmux passthrough (false)
%     'notes'          append the colour range and keys (true)
%   Any other name/value pair is a render option, see mio.render.

    [o, pairs] = mio.internal.render_split(varargin, {'cols', 80, 'rows', 24, ...
        'encoding', 'halfblock', 'color_depth', 'truecolor', 'format', 'ansi', ...
        'cell_aspect', 2, 'tmux', false, 'notes', true}, 'mio.render_text', true);
    bytes = meshioplusplus_mex('mio_render_text', mesh, pairs, mio.internal.render_text_opts(o));
    txt = mio.internal.render_utf8(bytes);
end
