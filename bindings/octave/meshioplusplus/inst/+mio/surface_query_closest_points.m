function r = surface_query_closest_points(q, points)
%SURFACE_QUERY_CLOSEST_POINTS  The nearest surface point of each query point.
%   R = mio.surface_query_closest_points(Q, POINTS) with POINTS (n x 3)
%   returns a struct: point (n x 3), distance (n x 1, unsigned) and cell
%   (n x 1, 1-based global input cell).

    r = meshioplusplus_mex('mio_surface_query_closest_points', q, double(points));
end
