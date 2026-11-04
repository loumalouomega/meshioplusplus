function q = surface_watertight_check(mesh)
%SURFACE_WATERTIGHT_CHECK  Whether a triangle surface is closed and consistently wound.
%   Q = mio.surface_watertight_check(MESH) returns a struct with
%   boundary_edges, non_manifold_edges, inconsistent_pairs,
%   degenerate_triangles and watertight (logical).

    q = meshioplusplus_mex('mio_surface_watertight_check', mesh);
end
