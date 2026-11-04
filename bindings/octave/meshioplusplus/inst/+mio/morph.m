function r = morph(mesh, array, varargin)
%MORPH  Boundary-driven mesh morphing by the point-data displacement ARRAY.
%   The constrained points (the boundary unless 'fix_boundary' is false, and
%   the nonzero entries of 'constrained_array') move by ARRAY; the rest by
%   its harmonic extension with Jacobian-based stiffening. Options:
%   'fix_boundary' (true), 'constrained_array' ([]), 'stiffening' (1),
%   'steps' (1) and 'allow_inverted' (false). Returns a struct: mesh,
%   num_constrained, num_free, num_inverted, max_displacement and converged.

    o = mio.internal.opts(varargin, {'fix_boundary', true, 'constrained_array', [], ...
        'stiffening', 1, 'steps', 1, 'allow_inverted', false}, 'mio.morph');
    r = meshioplusplus_mex('mio_morph', mesh, mio.internal.str(array, 'array'), ...
        logical(o.fix_boundary), mio.internal.str(o.constrained_array, 'constrained_array'), ...
        double(o.stiffening), double(o.steps), logical(o.allow_inverted));
end
