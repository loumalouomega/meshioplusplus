function f = render(mesh, varargin)
%RENDER  Render a mesh into an RGBA frame.
%   F = mio.render(MESH, 'name', value, ...) runs a deterministic software
%   rasterizer: no display, GPU or windowing library, and the same pixels on
%   every platform. Volume cells are drawn through their boundary skin. See
%   doc/tui.md. The options may also come as one struct with a field per
%   option. A [] value leaves that option at its default.
%
%   The options are those of the C mio_render_opts: 'width', 'height',
%   'supersample'; the camera 'view' ('iso', '+x' ... '-z'), 'azimuth',
%   'elevation', 'roll', 'projection' ('orthographic', 'perspective'),
%   'fov_deg', 'zoom', 'pan' ([x y]); 'shading' ('none', 'flat', 'smooth'),
%   'two_sided', 'ambient', 'light_dir' ([x y z]), 'split_angle'; 'edges'
%   ('none', 'all', 'feature'), 'feature_angle'; the colours 'edge_color',
%   'fill_color', 'line_color', 'background', 'nan_color', 'iso_color',
%   'vector_color', 'outline_color' (a 0xRRGGBBAA number from mio.rgba, or
%   [r g b] / [r g b a] in 0..255), 'point_radius'; colouring by 'color_by',
%   'component', 'cmap', 'vmin', 'vmax', 'colorbar', 'axes', 'scale_bar'; the
%   field options 'expr', 'reduce', 'clip' ([low high] percentiles, NaN for
%   either end), 'symmetric', 'scale' ('linear', 'log', 'symlog'),
%   'scale_threshold', 'categorical', 'color_regions', 'category_edges',
%   'isolines', 'iso_levels', 'vectors', 'vector_count', 'vector_length',
%   'warp', 'warp_scale', 'warp_outline', 'diagnostic' ('quality',
%   'inverted', 'degenerate', 'orientation', 'free_edges', 'edge_length'),
%   'quality_metric'; the cut-aways 'cutaway' (one or two planes of six
%   numbers: a point, then the normal of the side kept, as a vector of 6 or
%   12 numbers or a cell of 6-vectors) and 'cutaway_tint'; and the
%   streamlines 'streamlines' (a vector point array), 'stream_seeds',
%   'stream_length' and 'stream_color'. An unknown name is an error.
%
%   F is a struct:
%     width, height  frame size in pixels
%     rgba           height x width x 4 uint8, straight alpha. Row 1 is the
%                    TOP of the picture and column 1 its left edge, so
%                    image(F.rgba(:,:,1:3)) and imwrite(F.rgba(:,:,1:3),
%                    file, 'Alpha', F.rgba(:,:,4)) show it upright.
%     cell_ids       height x width double in the same orientation: the
%                    1-based global (block-major) cell drawn at each pixel,
%                    0 for none
%     range          [vmin vmax] of the mapped field, or [] when none
%     notes          n x 1 cellstr: the colour range, ticks and keys
%   mio.frame_array(F) gives the picture as doubles in [0, 1].

    [~, pairs] = mio.internal.render_split(varargin, {}, 'mio.render', false);
    f = meshioplusplus_mex('mio_render', mesh, pairs);
    notes = f.notes;
    for k = 1:numel(notes)
        notes{k} = mio.internal.render_utf8(notes{k});
    end
    f.notes = notes;
end
