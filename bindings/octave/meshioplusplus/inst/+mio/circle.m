function m = circle(varargin)
%CIRCLE  A closed circle of line cells.
%   M = mio.circle(...) takes these options, positionally in this order or
%   as name/value pairs (defaults in parentheses):
%     'radius'    radius (1)
%     'n'         number of cells (16)
%     'center'    centre of the shape ([0 0 0])
%     'order'     1; 2 quadratic on the exact geometry; 3 to 10 VTK_LAGRANGE_* (1)
%     'max_cells' refuse a mesh larger than this, by name (20000000)
%     'complete'  order 2: full-Lagrange quad9/hexahedron27/wedge18 (false)
%   Returns a mio.MeshHandle with named side regions and exact mid-edge
%   nodes for order 2. See doc/primitives.md for the point and cell order.

    o = mio.internal.opts(varargin, {'radius', 1, 'n', 16, 'center', [0 0 0], 'order', ...
        1, 'max_cells', 20000000, 'complete', false}, 'mio.circle');
    m = mio.internal.od_generate('circle', o.n, 'origin', o.center, 'radius', o.radius, ...
        'order', o.order, 'max_cells', o.max_cells, 'complete', o.complete);
end
