function r = laplacian_eigenmodes(mesh, varargin)
%LAPLACIAN_EIGENMODES  Lowest Laplace-Beltrami eigenpairs of a triangle surface.
%   The mesh carries point data 'laplacian:modes' (n x num_modes) and field
%   data 'laplacian:eigenvalues'. Options: 'num_modes' (10), 'tolerance'
%   (1e-8) and 'max_iterations' (0 takes 300). Returns a struct: mesh,
%   iterations, max_residual and converged.

    o = mio.internal.opts(varargin, {'num_modes', 10, 'tolerance', 1e-8, ...
        'max_iterations', 0}, 'mio.laplacian_eigenmodes');
    r = meshioplusplus_mex('mio_laplacian_eigenmodes', mesh, double(o.num_modes), ...
        double(o.tolerance), double(o.max_iterations));
end
