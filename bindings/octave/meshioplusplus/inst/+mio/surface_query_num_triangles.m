function n = surface_query_num_triangles(q)
%SURFACE_QUERY_NUM_TRIANGLES  The triangles of a surface query (quads fanned).

    n = meshioplusplus_mex('mio_surface_query_num_triangles', q);
end
