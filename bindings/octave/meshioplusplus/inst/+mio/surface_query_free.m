function surface_query_free(q)
%SURFACE_QUERY_FREE  Release a surface query handle now.
%   mio.surface_query_free(Q) frees the native query behind Q. Idempotent,
%   and a no-op for anything that is not a mio.SurfaceQueryHandle.

    meshioplusplus_mex('mio_surface_query_free', q);
end
