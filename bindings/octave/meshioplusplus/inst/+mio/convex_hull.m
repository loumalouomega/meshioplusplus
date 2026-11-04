function r = convex_hull(mesh, varargin)
%CONVEX_HULL  Convex hull of MESH's points.
%   In 3-D a closed outward triangle surface, in 2-D a counter-clockwise
%   line loop; each vertex carries 'hull:point_id' (its 0-based input
%   point). Options: 'dimension' (0 decides; 2 or 3). Returns a struct:
%   mesh, dimension, num_vertices, num_facets, area, volume and perimeter.

    o = mio.internal.opts(varargin, {'dimension', 0}, 'mio.convex_hull');
    r = meshioplusplus_mex('mio_convex_hull', mesh, double(o.dimension));
end
