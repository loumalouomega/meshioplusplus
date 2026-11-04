function m = plate_with_hole(varargin)
%PLATE_WITH_HOLE  A square plate with a central circular hole.
%   M = mio.plate_with_hole(...) takes these options, positionally in this order or
%   as name/value pairs (defaults in parentheses):
%     'radius'    hole radius (0.25)
%     'size'      plate side length (1)
%     'nt'        cells around the hole per 90 degree block (8)
%     'nr'        cells from the hole to the edge (4)
%     'center'    centre of the shape ([0 0 0])
%     'cell_type' 'quad', 'triangle', 'hexahedron', 'tetra' or 'wedge' as allowed; '' = default ('quad')
%     'order'     1; 2 quadratic on the exact geometry; 3 to 10 VTK_LAGRANGE_* (1)
%     'max_cells' refuse a mesh larger than this, by name (20000000)
%     'complete'  order 2: full-Lagrange quad9/hexahedron27/wedge18 (false)
%   Returns a mio.MeshHandle with named side regions and exact mid-edge
%   nodes for order 2. See doc/primitives.md for the point and cell order.

    o = mio.internal.opts(varargin, {'radius', 0.25, 'size', 1, 'nt', 8, 'nr', 4, ...
        'center', [0 0 0], 'cell_type', 'quad', 'order', 1, 'max_cells', 20000000, ...
        'complete', false}, 'mio.plate_with_hole');
    m = mio.internal.od_generate('plate_with_hole', [double(o.nt), double(o.nr)], ...
        'size', mio.internal.od_rep3(o.size), 'origin', o.center, 'radius', o.radius, ...
        'cell_type', o.cell_type, 'order', o.order, 'max_cells', o.max_cells, ...
        'complete', o.complete);
end
