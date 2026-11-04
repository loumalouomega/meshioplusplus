function m = torus(varargin)
%TORUS  A torus surface.
%   M = mio.torus(...) takes these options, positionally in this order or
%   as name/value pairs (defaults in parentheses):
%     'major_radius' major radius (1)
%     'minor_radius' minor radius (0.25)
%     'nu'           cells around the axis (32)
%     'nv'           cells around the tube (16)
%     'center'       centre of the shape ([0 0 0])
%     'cell_type'    'quad', 'triangle', 'hexahedron', 'tetra' or 'wedge' as allowed; '' = default ('quad')
%     'order'        1; 2 quadratic on the exact geometry; 3 to 10 VTK_LAGRANGE_* (1)
%     'max_cells'    refuse a mesh larger than this, by name (20000000)
%     'complete'     order 2: full-Lagrange quad9/hexahedron27/wedge18 (false)
%   Returns a mio.MeshHandle with named side regions and exact mid-edge
%   nodes for order 2. See doc/primitives.md for the point and cell order.

    o = mio.internal.opts(varargin, {'major_radius', 1, 'minor_radius', 0.25, 'nu', 32, ...
        'nv', 16, 'center', [0 0 0], 'cell_type', 'quad', 'order', 1, 'max_cells', ...
        20000000, 'complete', false}, 'mio.torus');
    m = mio.internal.od_generate('torus', [double(o.nu), double(o.nv)], 'origin', ...
        o.center, 'radius', o.major_radius, 'radius2', o.minor_radius, 'cell_type', ...
        o.cell_type, 'order', o.order, 'max_cells', o.max_cells, 'complete', o.complete);
end
