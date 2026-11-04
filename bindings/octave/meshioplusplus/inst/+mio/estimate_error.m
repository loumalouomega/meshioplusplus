function r = estimate_error(mesh, array, varargin)
%ESTIMATE_ERROR  Recovery-based (Zienkiewicz-Zhu) error indicator of ARRAY.
%   Options: 'method' ('zz'), 'marking' ('none', 'absolute', 'fraction' or
%   'dorfler'), 'marking_value' (0), 'output' ('' uses 'error:zz'),
%   'marked' ('' uses 'error:marked') and 'overwrite' (false). Returns a
%   struct: mesh, global_error, num_skipped and num_marked.

    o = mio.internal.opts(varargin, {'method', 'zz', 'marking', 'none', 'marking_value', 0, ...
        'output', '', 'marked', '', 'overwrite', false}, 'mio.estimate_error');
    r = meshioplusplus_mex('mio_estimate_error', mesh, mio.internal.str(array, 'array'), ...
        mio.internal.str(o.method, 'method'), mio.internal.str(o.marking, 'marking'), ...
        double(o.marking_value), mio.internal.str(o.output, 'output'), ...
        mio.internal.str(o.marked, 'marked'), logical(o.overwrite));
end
