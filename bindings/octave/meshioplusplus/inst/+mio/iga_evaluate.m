function ev = iga_evaluate(info, patch, params, varargin)
%IGA_EVALUATE  Evaluate a NURBS patch at parametric coordinates.
%   EV = mio.iga_evaluate(INFO, PATCH, PARAMS, ...) evaluates patch PATCH
%   (1-based; the face of a B-rep) of the model in the mio.FormatInfo INFO
%   at the rows of the (N x d) matrix PARAMS, one row per point, each inside
%   the knot range (d the patch's parametric dimension, 2 for a face).
%   Options (positional in this order, or name/value):
%     'derivatives'   0, or 1 for the first derivatives and normals (0)
%     'trim_samples'  B-rep: points per knot span a trim is sampled at (16)
%   EV is a struct with
%     points       N x space_dimension
%     fields       cell, one N x components matrix per NURBS field
%     derivatives  N x d x space_dimension (d x / d u_k), or []
%     normals      N x 3 unit normals of a surface in 3-D, or []
%     inside       N x 1 logical for a B-rep face, or []
%   derivatives and normals are [] unless 'derivatives' is 1.

    o = mio.internal.opts(varargin, {'derivatives', 0, 'trim_samples', 16}, 'mio.iga_evaluate');
    ev = meshioplusplus_mex('mio_iga_evaluate', info, double(patch), double(params), ...
        double(o.derivatives), double(o.trim_samples));
end
