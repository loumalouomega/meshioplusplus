function m = reference_triangle(varargin)
%REFERENCE_TRIANGLE  The reference right triangle, subdivided.
%   M = mio.reference_triangle(...) takes these options, positionally in this order or
%   as name/value pairs (defaults in parentheses):
%     'n'         divisions per side (4)
%     'size'      leg lengths ([1 1])
%     'origin'    right-angle corner ([0 0 0])
%     'order'     1; 2 quadratic on the exact geometry; 3 to 10 VTK_LAGRANGE_* (1)
%     'max_cells' refuse a mesh larger than this, by name (20000000)
%     'complete'  order 2: full-Lagrange quad9/hexahedron27/wedge18 (false)
%   Returns a mio.MeshHandle with named side regions and exact mid-edge
%   nodes for order 2. See doc/primitives.md for the point and cell order.

    o = mio.internal.opts(varargin, {'n', 4, 'size', [1 1], 'origin', [0 0 0], 'order', ...
        1, 'max_cells', 20000000, 'complete', false}, 'mio.reference_triangle');
    m = mio.internal.od_generate('reference_triangle', o.n, 'size', ...
        mio.internal.od_pad3(o.size), 'origin', o.origin, 'order', o.order, ...
        'max_cells', o.max_cells, 'complete', o.complete);
end
