function m = annulus(varargin)
%ANNULUS  An annulus, or an annular sector.
%   M = mio.annulus(...) takes these options, positionally in this order or
%   as name/value pairs (defaults in parentheses):
%     'inner_radius' inner radius (0.5)
%     'outer_radius' outer radius (1)
%     'ntheta'       cells around (16)
%     'nr'           cells across the gap (4)
%     'angle'        opening angle, radians; a full ring welds its seam (2 * pi)
%     'center'       centre of the shape ([0 0 0])
%     'cell_type'    'quad', 'triangle', 'hexahedron', 'tetra' or 'wedge' as allowed; '' = default ('quad')
%     'order'        1; 2 quadratic on the exact geometry; 3 to 10 VTK_LAGRANGE_* (1)
%     'max_cells'    refuse a mesh larger than this, by name (20000000)
%     'complete'     order 2: full-Lagrange quad9/hexahedron27/wedge18 (false)
%   Returns a mio.MeshHandle with named side regions and exact mid-edge
%   nodes for order 2. See doc/primitives.md for the point and cell order.

    o = mio.internal.opts(varargin, {'inner_radius', 0.5, 'outer_radius', 1, 'ntheta', ...
        16, 'nr', 4, 'angle', 2 * pi, 'center', [0 0 0], 'cell_type', 'quad', 'order', ...
        1, 'max_cells', 20000000, 'complete', false}, 'mio.annulus');
    m = mio.internal.od_generate('annulus', [double(o.ntheta), double(o.nr)], 'origin', ...
        o.center, 'radius', o.outer_radius, 'radius2', o.inner_radius, 'angle', ...
        o.angle, 'cell_type', o.cell_type, 'order', o.order, 'max_cells', o.max_cells, ...
        'complete', o.complete);
end
