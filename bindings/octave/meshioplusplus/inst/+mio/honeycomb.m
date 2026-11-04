function m = honeycomb(varargin)
%HONEYCOMB  A planar honeycomb of hexagonal polygons.
%   M = mio.honeycomb(...) takes these options, positionally in this order or
%   as name/value pairs (defaults in parentheses):
%     'dims'      hexagons [nx ny] ([4 4])
%     'edge'      side length of a hexagon (1)
%     'origin'    lo corner ([0 0 0])
%     'max_cells' refuse a mesh larger than this, by name (20000000)
%   Returns a mio.MeshHandle with named side regions and exact mid-edge
%   nodes for order 2. See doc/primitives.md for the point and cell order.

    o = mio.internal.opts(varargin, {'dims', [4 4], 'edge', 1, 'origin', [0 0 0], ...
        'max_cells', 20000000}, 'mio.honeycomb');
    m = mio.internal.od_generate('honeycomb', o.dims, 'size', ...
        mio.internal.od_rep3(o.edge), 'origin', o.origin, 'max_cells', o.max_cells);
end
