function r = iga_remove_knots(info, varargin)
%IGA_REMOVE_KNOTS  Remove knots from a NURBS patch (The NURBS Book, A5.8).
%   R = mio.iga_remove_knots(INFO, ...) removes interior knots from one
%   direction of one patch of the NURBS model in the mio.FormatInfo INFO while
%   the geometry and the NURBS fields move by at most the tolerance. Options
%   (positional in this order, or name/value):
%     'direction'            the parametric direction, 1-based (1)
%     'knots'                knot values to remove, one copy per entry; [] is
%                            every copy of every interior knot ([])
%     'tolerance'            the largest deviation; 0 removes only what is exact
%                            up to rounding, which undoes a knot insertion (0)
%     'patch'                the patch, 1-based (1)
%     'sample_subdivisions'  cells per knot span of the returned sampling (0)
%   R is a struct: 'mesh', 'info' (the reduced model) and 'report', a struct
%   with 'knots', 'requested' and 'removed' (per distinct knot value, as
%   column vectors) and 'num_removed'. A B-rep, and a model whose patch
%   topology joins several patches, are refused.

    o = mio.internal.opts(varargin, {'direction', 1, 'knots', [], 'tolerance', 0, ...
        'patch', 1, 'sample_subdivisions', 0}, 'mio.iga_remove_knots');
    r = meshioplusplus_mex('mio_iga_remove_knots', info, double(o.patch), double(o.direction), ...
        double(o.knots(:)), double(o.tolerance), double(o.sample_subdivisions));
end
