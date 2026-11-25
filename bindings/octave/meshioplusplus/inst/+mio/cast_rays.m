function r = cast_rays(mesh, origins, directions, varargin)
%CAST_RAYS  The nearest surface hit of a batch of rays.
%   R = mio.cast_rays(MESH, ORIGINS, DIRECTIONS) shoots every ray
%   origin + t * direction (the direction is normalised, so t is a length) at
%   the surface cells of MESH (triangles, quadrilaterals and polygons) and
%   reports the nearest hit with 0 < t <= max_distance. ORIGINS is n x 3, one
%   ray per row; DIRECTIONS is 1 x 3 (one direction for every ray) or n x 3.
%   Options: 'max_distance' (0: no limit) and 'region' ('': every surface cell).
%   Returns a struct: hit (n x 1 logical), distance (n x 1, NaN on a miss),
%   point and normal (n x 3, unit normal on the triangle's winding side, NaN on
%   a miss), cell (n x 1, the 1-based input cell hit, 0 on a miss) and
%   num_hits.

    o = mio.internal.opts(varargin, {'max_distance', 0, 'region', ''}, 'mio.cast_rays');
    r = meshioplusplus_mex('mio_cast_rays', mesh, double(origins), double(directions), ...
        double(o.max_distance), mio.internal.str(o.region, 'region'));
end
