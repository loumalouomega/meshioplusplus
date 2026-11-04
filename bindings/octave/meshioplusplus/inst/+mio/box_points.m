function m = box_points(varargin)
%BOX_POINTS  Pseudo-random points in a box (vertex cells).
%   M = mio.box_points(...) takes these options, positionally in this order or
%   as name/value pairs (defaults in parentheses):
%     'n'         number of points (100)
%     'size'      side lengths ([1 1 1])
%     'origin'    lo corner ([0 0 0])
%     'seed'      the splitmix64 seed (0)
%     'max_cells' refuse a mesh larger than this, by name (20000000)
%   Returns a mio.MeshHandle with named side regions and exact mid-edge
%   nodes for order 2. See doc/primitives.md for the point and cell order.

    o = mio.internal.opts(varargin, {'n', 100, 'size', [1 1 1], 'origin', [0 0 0], ...
        'seed', 0, 'max_cells', 20000000}, 'mio.box_points');
    m = mio.internal.od_generate('box_points', o.n, 'size', o.size, 'origin', o.origin, ...
        'seed', o.seed, 'max_cells', o.max_cells);
end
