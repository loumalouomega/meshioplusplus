function r = extrude(mesh, varargin)
%EXTRUDE  Sweep a surface (or a curve) along a vector.
%   R = mio.extrude(MESH, ...) sweeps the triangles and quads of MESH into
%   wedges and hexahedra (or, with no 2-D cells, its lines into quads) along
%   'direction', the whole displacement, in 'num_layers' layers spaced
%   uniformly, by the geometric 'ratio', or at the explicit 'offsets'
%   (fractions, increasing, ending at 1) (roadmap 3.1.2). The result holds
%   the swept cells only, with Side regions 'bottom', 'top' and 'side'.
%   Options:
%     'direction'     the whole displacement [x y z] ([0 0 1])
%     'num_layers'    layers (1)
%     'ratio'         each layer's thickness over the one below's (1)
%     'offsets'       explicit level fractions, or empty (zeros(0, 1))
%     'record_layers' attach 'extrude:layer' cell data (false)
%     'max_cells'     refuse more profile cells times layers (20000000)
%   R is a struct: mesh, num_cells, num_hexahedra, num_wedges, num_pyramids,
%   num_tetrahedra, num_quads, num_triangles, num_flipped, num_on_axis,
%   min_scaled_jacobian, full and warnings (n x 1 cellstr).

    o = mio.internal.opts(varargin, {'direction', [0 0 1], 'num_layers', 1, 'ratio', 1, ...
        'offsets', zeros(0, 1), 'record_layers', false, 'max_cells', 20000000}, 'mio.extrude');
    r = meshioplusplus_mex('mio_extrude', mesh, double(o.direction), double(o.num_layers), ...
        double(o.ratio), double(o.offsets), logical(o.record_layers), double(o.max_cells));
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
