function r = iga_closest_point(info, points, varargin)
%IGA_CLOSEST_POINT  The point of a NURBS patch closest to each query point.
%   R = mio.iga_closest_point(INFO, POINTS) finds, for each row of the (N x
%   space_dim) matrix POINTS, the point of patch 1 of the NURBS model in the
%   mio.FormatInfo INFO closest to it (a curve or a surface, rational or not).
%   R = mio.iga_closest_point(INFO, POINTS, PATCH) uses patch PATCH (1-based).
%   A point beyond the patch's edge lands on the edge; a volume is refused.
%   R is a struct with
%     params     N x d parameters of the closest points
%     points     N x space_dimension closest points
%     distances  N x 1 distances to them

    if isempty(varargin)
        patch = 1;
    else
        patch = varargin{1};
    end
    r = meshioplusplus_mex('mio_iga_closest_point', info, double(patch), double(points));
end
