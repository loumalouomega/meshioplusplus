function r = compute_curvature(mesh, varargin)
%COMPUTE_CURVATURE  Per-vertex mean and Gaussian curvature of a surface.
%   Writes 'curvature:mean' and 'curvature:gaussian' point data, optionally
%   'curvature:area' and 'curvature:principal' (n x 2, k1 >= k2). Options:
%   'mean' (true), 'gaussian' (true), 'dual_area' ('mixed-voronoi' or
%   'barycentric'), 'include_boundary' (false), 'record_area' (false),
%   'record_principal' (false) and 'region' ('' for every surface cell).
%   Returns a struct: mesh, quality (a struct of boundary_edges,
%   non_manifold_edges, inconsistent_pairs, degenerate_triangles,
%   watertight), num_boundary, num_isolated, num_degenerate and
%   total_angle_defect (4*pi on any closed sphere-like surface).

    o = mio.internal.opts(varargin, {'mean', true, 'gaussian', true, ...
        'dual_area', 'mixed-voronoi', 'include_boundary', false, 'record_area', false, ...
        'record_principal', false, 'region', ''}, 'mio.compute_curvature');
    r = meshioplusplus_mex('mio_compute_curvature', mesh, logical(o.mean), logical(o.gaussian), ...
        mio.internal.str(o.dual_area, 'dual_area'), logical(o.include_boundary), ...
        logical(o.record_area), logical(o.record_principal), mio.internal.str(o.region, 'region'));
end
