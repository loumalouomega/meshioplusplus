function r = revolve(mesh, varargin)
%REVOLVE  Sweep a surface (or a curve) around an axis.
%   R = mio.revolve(MESH, ...) turns the triangles and quads of MESH (or its
%   lines) around the axis through 'axis_point' along 'axis' by 'angle'
%   degrees in 'num_layers' steps (roadmap 3.1.3); points within 'tolerance'
%   of the axis stay put and the cells around them collapse. The result holds
%   the swept cells only, with Side regions 'bottom', 'top' ('start', 'end'
%   for a partial turn) and 'side'. Options:
%     'axis'          the axis direction ([0 1 0])
%     'axis_point'    a point of the axis ([0 0 0])
%     'angle'         degrees, right-hand rule; 360 closes the turn (360)
%     'num_layers'    steps of the turn (16)
%     'tolerance'     on the axis: a fraction of the bounding-box diagonal (1e-9)
%     'record_layers' attach 'revolve:layer' cell data (false)
%     'max_cells'     refuse more profile cells times layers (20000000)
%   R is the struct mio.extrude returns.

    o = mio.internal.opts(varargin, {'axis', [0 1 0], 'axis_point', [0 0 0], 'angle', 360, ...
        'num_layers', 16, 'tolerance', 1e-9, 'record_layers', false, ...
        'max_cells', 20000000}, 'mio.revolve');
    r = meshioplusplus_mex('mio_revolve', mesh, double(o.axis), double(o.axis_point), ...
        double(o.angle), double(o.num_layers), double(o.tolerance), logical(o.record_layers), ...
        double(o.max_cells));
    r.warnings = split_warnings(r.warnings);
end

function w = split_warnings(text)
    if isempty(text)
        w = cell(0, 1);
    else
        w = strsplit(text, char(10));
        w = w(:);
    end
end
