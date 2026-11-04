function m = spherical_shell(varargin)
%SPHERICAL_SHELL  A solid spherical shell.
%   M = mio.spherical_shell(...) takes these options, positionally in this order or
%   as name/value pairs (defaults in parentheses):
%     'inner_radius' inner radius (0.5)
%     'outer_radius' outer radius (1)
%     'n'            cells per cube face side (4)
%     'nr'           layers (2)
%     'ratio'        layer thickness ratio (1)
%     'center'       centre of the shape ([0 0 0])
%     'order'        1; 2 quadratic on the exact geometry; 3 to 10 VTK_LAGRANGE_* (1)
%     'max_cells'    refuse a mesh larger than this, by name (20000000)
%     'complete'     order 2: full-Lagrange quad9/hexahedron27/wedge18 (false)
%   Returns a mio.MeshHandle with named side regions and exact mid-edge
%   nodes for order 2. See doc/primitives.md for the point and cell order.

    o = mio.internal.opts(varargin, {'inner_radius', 0.5, 'outer_radius', 1, 'n', 4, ...
        'nr', 2, 'ratio', 1, 'center', [0 0 0], 'order', 1, 'max_cells', 20000000, ...
        'complete', false}, 'mio.spherical_shell');
    m = mio.internal.od_generate('spherical_shell', [double(o.n), double(o.nr)], ...
        'origin', o.center, 'radius', o.outer_radius, 'radius2', o.inner_radius, ...
        'ratio', o.ratio, 'order', o.order, 'max_cells', o.max_cells, 'complete', ...
        o.complete);
end
