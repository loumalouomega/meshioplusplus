function png = render_png(mesh, varargin)
%RENDER_PNG  Render a mesh to PNG bytes.
%   PNG = mio.render_png(MESH, ...) returns an RGBA PNG file as an n x 1
%   uint8 column; fwrite(fid, PNG) saves it.
%   Options, positional or as a name/value pair:
%     'compress'   0 (default) stores the pixels uncompressed (no zlib, the
%                  same bytes everywhere); 1-9 compress through zlib where
%                  the build has it
%   Any other name/value pair is a render option, see mio.render.

    [o, pairs] = mio.internal.render_split(varargin, {'compress', 0}, 'mio.render_png', true);
    png = meshioplusplus_mex('mio_render_png', mesh, pairs, double(o.compress));
end
