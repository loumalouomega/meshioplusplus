function m = sphere_points(varargin)
%SPHERE_POINTS  Points spread evenly over a sphere (vertex cells).
%   M = mio.sphere_points(...) takes these options, positionally in this order or
%   as name/value pairs (defaults in parentheses):
%     'n'         number of points (100)
%     'radius'    radius (1)
%     'center'    centre of the shape ([0 0 0])
%     'max_cells' refuse a mesh larger than this, by name (20000000)
%   Returns a mio.MeshHandle with named side regions and exact mid-edge
%   nodes for order 2. See doc/primitives.md for the point and cell order.

    o = mio.internal.opts(varargin, {'n', 100, 'radius', 1, 'center', [0 0 0], ...
        'max_cells', 20000000}, 'mio.sphere_points');
    m = mio.internal.od_generate('sphere_points', o.n, 'origin', o.center, 'radius', ...
        o.radius, 'max_cells', o.max_cells);
end
