function m = disk(varargin)
%DISK  A disk of quads or triangles.
%   M = mio.disk(...) takes these options, positionally in this order or
%   as name/value pairs (defaults in parentheses):
%     'radius'    radius (1)
%     'n'         cells across the diameter (quads) or number of rings (triangles) (8)
%     'center'    centre of the shape ([0 0 0])
%     'cell_type' 'quad', 'triangle', 'hexahedron', 'tetra' or 'wedge' as allowed; '' = default ('quad')
%     'order'     1; 2 quadratic on the exact geometry; 3 to 10 VTK_LAGRANGE_* (1)
%     'max_cells' refuse a mesh larger than this, by name (20000000)
%     'complete'  order 2: full-Lagrange quad9/hexahedron27/wedge18 (false)
%   Returns a mio.MeshHandle with named side regions and exact mid-edge
%   nodes for order 2. See doc/primitives.md for the point and cell order.

    o = mio.internal.opts(varargin, {'radius', 1, 'n', 8, 'center', [0 0 0], ...
        'cell_type', 'quad', 'order', 1, 'max_cells', 20000000, 'complete', false}, ...
        'mio.disk');
    m = mio.internal.od_generate('disk', o.n, 'origin', o.center, 'radius', o.radius, ...
        'cell_type', o.cell_type, 'order', o.order, 'max_cells', o.max_cells, ...
        'complete', o.complete);
end
