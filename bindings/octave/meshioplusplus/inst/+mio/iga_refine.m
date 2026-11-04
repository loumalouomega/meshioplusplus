function r = iga_refine(info, varargin)
%IGA_REFINE  Refine a NURBS model exactly (the isogeometric layer).
%   R = mio.iga_refine(INFO, ...) refines every patch (every face surface of
%   a B-rep) of the NURBS model in the mio.FormatInfo INFO (an MFEM NURBS
%   mesh, Kratos CAD JSON, STEP or IGES, from mio.read_with_info): degree
%   elevation by 'elevate', then knot insertion splitting every nonzero knot
%   span in 'subdivisions' parts. The geometry, the fields and a trimmed
%   face's trimming loops do not change.
%   Options (positional in this order, or name/value):
%     'subdivisions'         one integer, or one per parametric direction (1)
%     'elevate'              one integer, or one per direction (0)
%     'sample_subdivisions'  cells per knot span of the sampling; 0 is
%                            automatic (0)
%     'trim_samples'         B-rep: points per knot span a trimming curve
%                            is sampled at (16)
%     'edges'                B-rep: sample the edge curves (true)
%     'trim'                 B-rep: keep only cells inside the trims (true)
%   R is a struct: 'mesh' (the refined model sampled) and 'info' (the
%   refined model, for mio.write_with_info).

    o = mio.internal.opts(varargin, {'subdivisions', 1, 'elevate', 0, ...
        'sample_subdivisions', 0, 'trim_samples', 16, 'edges', true, 'trim', true}, ...
        'mio.iga_refine');
    r = meshioplusplus_mex('mio_iga_refine', info, double(o.subdivisions), double(o.elevate), ...
        double(o.sample_subdivisions), double(o.trim_samples), logical(o.edges), logical(o.trim));
end
