function m = line(varargin)
%LINE  A straight polyline of line cells.
%   M = mio.line(...) takes these options, positionally in this order or
%   as name/value pairs (defaults in parentheses):
%     'n'         number of cells (8)
%     'start'     first end point ([0 0 0])
%     'stop'      second end point (R's end, a reserved word here) ([1 0 0])
%     'order'     1; 2 quadratic on the exact geometry; 3 to 10 VTK_LAGRANGE_* (1)
%     'max_cells' refuse a mesh larger than this, by name (20000000)
%     'complete'  order 2: full-Lagrange quad9/hexahedron27/wedge18 (false)
%   Returns a mio.MeshHandle with named side regions and exact mid-edge
%   nodes for order 2. See doc/primitives.md for the point and cell order.

    o = mio.internal.opts(varargin, {'n', 8, 'start', [0 0 0], 'stop', [1 0 0], ...
        'order', 1, 'max_cells', 20000000, 'complete', false}, 'mio.line');
    m = mio.internal.od_generate('line', o.n, 'size', o.stop, 'origin', o.start, ...
        'order', o.order, 'max_cells', o.max_cells, 'complete', o.complete);
end
