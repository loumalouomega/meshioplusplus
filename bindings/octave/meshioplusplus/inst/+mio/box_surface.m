function m = box_surface(varargin)
%BOX_SURFACE  The closed surface of a box.
%   M = mio.box_surface(...) takes these options, positionally in this order or
%   as name/value pairs (defaults in parentheses):
%     'dims'      cell counts [nx ny nz] ([1 1 1])
%     'size'      side lengths ([1 1 1])
%     'origin'    lo corner ([0 0 0])
%     'cell_type' 'quad', 'triangle', 'hexahedron', 'tetra' or 'wedge' as allowed; '' = default ('quad')
%     'order'     1; 2 quadratic on the exact geometry; 3 to 10 VTK_LAGRANGE_* (1)
%     'max_cells' refuse a mesh larger than this, by name (20000000)
%     'complete'  order 2: full-Lagrange quad9/hexahedron27/wedge18 (false)
%   Returns a mio.MeshHandle with named side regions and exact mid-edge
%   nodes for order 2. See doc/primitives.md for the point and cell order.

    o = mio.internal.opts(varargin, {'dims', [1 1 1], 'size', [1 1 1], 'origin', ...
        [0 0 0], 'cell_type', 'quad', 'order', 1, 'max_cells', 20000000, 'complete', ...
        false}, 'mio.box_surface');
    m = mio.internal.od_generate('box_surface', o.dims, 'size', o.size, 'origin', ...
        o.origin, 'cell_type', o.cell_type, 'order', o.order, 'max_cells', o.max_cells, ...
        'complete', o.complete);
end
