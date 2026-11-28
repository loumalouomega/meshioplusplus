function r = surface_query_cast_rays(q, origins, directions, varargin)
%SURFACE_QUERY_CAST_RAYS  mio.cast_rays against a surface query's surface.
%   R = mio.surface_query_cast_rays(Q, ORIGINS, DIRECTIONS, ...) with ORIGINS
%   (n x 3), DIRECTIONS (1 x 3 or n x 3) and 'max_distance' (0, no limit)
%   returns the struct mio.cast_rays does, bit-identically.

    o = mio.internal.opts(varargin, {'max_distance', 0}, 'mio.surface_query_cast_rays');
    r = meshioplusplus_mex('mio_surface_query_cast_rays', q, double(origins), ...
        double(directions), double(o.max_distance));
end
