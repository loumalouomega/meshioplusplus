function m = lshape(varargin)
%LSHAPE  An L-shaped planar domain.
%   M = mio.lshape(...) takes these options, positionally in this order or
%   as name/value pairs (defaults in parentheses):
%     'n'         cells per unit block side (4)
%     'size'      overall side length (2)
%     'center'    centre of the shape ([0 0 0])
%     'cell_type' 'quad', 'triangle', 'hexahedron', 'tetra' or 'wedge' as allowed; '' = default ('quad')
%     'order'     1; 2 quadratic on the exact geometry; 3 to 10 VTK_LAGRANGE_* (1)
%     'max_cells' refuse a mesh larger than this, by name (20000000)
%     'complete'  order 2: full-Lagrange quad9/hexahedron27/wedge18 (false)
%   Returns a mio.MeshHandle with named side regions and exact mid-edge
%   nodes for order 2. See doc/primitives.md for the point and cell order.

    o = mio.internal.opts(varargin, {'n', 4, 'size', 2, 'center', [0 0 0], 'cell_type', ...
        'quad', 'order', 1, 'max_cells', 20000000, 'complete', false}, 'mio.lshape');
    m = mio.internal.od_generate('lshape', o.n, 'size', mio.internal.od_rep3(o.size), ...
        'origin', o.center, 'cell_type', o.cell_type, 'order', o.order, 'max_cells', ...
        o.max_cells, 'complete', o.complete);
end
