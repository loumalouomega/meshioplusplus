function write_gltf(mesh, path, varargin)
%WRITE_GLTF  Export a coloured glTF surface (.glb or .gltf).
%   Options (R defaults): 'color_by' ([]), 'cmap' ('viridis'), 'component'
%   (1-based, [] for magnitude), 'vmin'/'vmax' ([] for automatic),
%   'split_angle' (30), 'up_axis' ('auto','z','y','x'), 'scale' (1),
%   'container' ('auto','glb','gltf'), 'normal_weight' ('angle','area'),
%   'normals', 'fields', 'recenter', 'by_region', 'unlit' (all true),
%   'nan_color' ('#808080').

    o = mio.internal.opts(varargin, {'color_by', [], 'cmap', 'viridis', 'component', [], ...
        'vmin', [], 'vmax', [], 'split_angle', 30, 'up_axis', 'auto', 'scale', 1, ...
        'container', 'auto', 'normal_weight', 'angle', 'normals', true, 'fields', true, ...
        'recenter', true, 'by_region', true, 'unlit', true, 'nan_color', '#808080'}, ...
        'mio.write_gltf');
    if isempty(o.component)
        component = 1;
    else
        component = double(o.component);
    end
    if isempty(o.vmin)
        vmin = 0;
    else
        vmin = double(o.vmin);
    end
    if isempty(o.vmax)
        vmax = 0;
    else
        vmax = double(o.vmax);
    end
    options = {mio.internal.choice(o.container, {'auto', 'glb', 'gltf'}, 'container'), ...
        mio.internal.choice(o.up_axis, {'auto', 'z', 'y', 'x'}, 'up_axis'), ...
        mio.internal.choice(o.normal_weight, {'angle', 'area'}, 'normal_weight'), ...
        logical(o.normals), logical(o.fields), logical(o.recenter), logical(o.by_region), ...
        logical(o.unlit), component, ~isempty(o.component), ~isempty(o.vmin), ...
        ~isempty(o.vmax), double(o.split_angle), double(o.scale), vmin, vmax, ...
        mio.internal.str(o.color_by, 'color_by'), mio.internal.str(o.cmap, 'cmap'), ...
        mio.internal.str(o.nan_color, 'nan_color')};
    meshioplusplus_mex('mio_write_gltf', mesh, mio.internal.str(path, 'path'), options);
end
