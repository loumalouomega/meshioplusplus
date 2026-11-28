function q = surface_query(mesh, varargin)
%SURFACE_QUERY  A reusable surface query handle (roadmap 9.3.2).
%   Q = mio.surface_query(MESH, ...) builds the triangle soup and hierarchies
%   of MESH's triangle, quad and polygon cells (or of cell region 'region')
%   once; mio.surface_query_cast_rays and mio.surface_query_closest_points
%   then query them. The handle copies the surface; mio.surface_query_free
%   releases it (or it is released when cleared). See doc/ray_queries.md.

    o = mio.internal.opts(varargin, {'region', []}, 'mio.surface_query');
    q = meshioplusplus_mex('mio_surface_query', mesh, mio.internal.str(o.region, 'region'));
end
