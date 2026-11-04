function r = sobolev_deform(mesh, array, length_scale, varargin)
%SOBOLEV_DEFORM  Move MESH's points by a Helmholtz-filtered displacement.
%   Solves (M + l^2 K) u = M d for the (n x dim) point-data displacement
%   ARRAY with l = LENGTH_SCALE. Options: 'fixed_points_array' ('': a
%   point-data array whose nonzero entries pin their point), 'fix_boundary'
%   (false), 'record_filtered' (false), 'max_iterations' (128) and
%   'tolerance' (1e-10). Returns a struct: mesh, num_iterations, residual,
%   converged, num_fixed, num_isolated and max_displacement.

    o = mio.internal.opts(varargin, {'fixed_points_array', '', 'fix_boundary', false, ...
        'record_filtered', false, 'max_iterations', 128, 'tolerance', 1e-10}, ...
        'mio.sobolev_deform');
    r = meshioplusplus_mex('mio_sobolev_deform', mesh, mio.internal.str(array, 'array'), ...
        double(length_scale), mio.internal.str(o.fixed_points_array, 'fixed_points_array'), ...
        logical(o.fix_boundary), logical(o.record_filtered), double(o.max_iterations), ...
        double(o.tolerance));
end
