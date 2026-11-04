function r = iga_sample(info, varargin)
%IGA_SAMPLE  Sample a NURBS model into Lagrange cells.
%   R = mio.iga_sample(INFO, ...) samples the NURBS model in the
%   mio.FormatInfo INFO. Options (positional in this order, or name/value):
%     'sample_subdivisions'  cells per knot span and direction; 0 is
%                            automatic (0)
%     'trim_samples'         B-rep: points per knot span of a trim (16)
%     'edges'                B-rep: sample the edge curves (true)
%     'trim'                 B-rep: keep only cells inside the trims (true)
%   R is a struct: 'mesh' and 'info' (the model with the sampling's
%   cell-to-patch maps).

    o = mio.internal.opts(varargin, {'sample_subdivisions', 0, 'trim_samples', 16, ...
        'edges', true, 'trim', true}, 'mio.iga_sample');
    r = meshioplusplus_mex('mio_iga_sample', info, double(o.sample_subdivisions), ...
        double(o.trim_samples), logical(o.edges), logical(o.trim));
end
