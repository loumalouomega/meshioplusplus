function r = colorize(mesh, field, varargin)
%COLORIZE  Bake a field into red, green and blue colours.
%   R = mio.colorize(MESH, FIELD, ...) maps FIELD (point or cell data) through
%   a colormap into the UInt8 data 'red', 'green' and 'blue' on the field's
%   location (roadmap 7.4.4). Options:
%     'cmap'      a colormap name ('viridis')
%     'vmin'      the range minimum, [] automatic ([])
%     'vmax'      the range maximum, [] automatic ([])
%     'symmetric' with the automatic range, +-max|v| (false)
%     'steps'     bands at fixed steps, 0 continuous (0)
%     'nan_color' '#rrggbb' for a non-finite value ('#808080')
%     'component' colour by this 1-based component, [] the magnitude ([])
%   R is a struct: mesh, vmin, vmax and location ('point' or 'cell').

    o = mio.internal.opts(varargin, {'cmap', 'viridis', 'vmin', [], 'vmax', [], ...
        'symmetric', false, 'steps', 0, 'nan_color', '#808080', 'component', []}, ...
        'mio.colorize');
    if isempty(o.vmin) ~= isempty(o.vmax)
        error('mio:colorize', 'mio.colorize: give both vmin and vmax, or neither');
    end
    lo = 0;
    hi = 0;
    if ~isempty(o.vmin)
        lo = double(o.vmin);
        hi = double(o.vmax);
    end
    comp = -1;
    if ~isempty(o.component)
        comp = double(o.component) - 1;
    end
    r = meshioplusplus_mex('mio_colorize', mesh, mio.internal.str(field, 'field'), ...
        mio.internal.str(o.cmap, 'cmap'), lo, hi, double(logical(o.symmetric)), ...
        double(o.steps), mio.internal.str(o.nan_color, 'nan_color'), comp);
end
