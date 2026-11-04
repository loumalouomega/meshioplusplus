function m = blend_steps(a, b, w, varargin)
%BLEND_STEPS  Blend two steps of one topology.
%   M = mio.blend_steps(A, B, W, ...) is A with every floating-point data
%   array replaced by (1 - W) A + W B; integer data comes from the nearer
%   step. A and B have the same point count and blocks. Option
%   'blend_points' (blend the point coordinates too; false). See
%   doc/sequences.md. M is a new mesh.

    o = mio.internal.opts(varargin, {'blend_points', false}, 'mio.blend_steps');
    m = meshioplusplus_mex('mio_blend_steps', a, b, double(w), logical(o.blend_points));
end
