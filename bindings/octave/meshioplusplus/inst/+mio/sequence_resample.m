function sequence_resample(seq, out_path, times, varargin)
%SEQUENCE_RESAMPLE  Resample a sequence onto new times.
%   mio.sequence_resample(SEQ, OUT_PATH, TIMES, ...) writes one file per
%   target time in TIMES when OUT_PATH carries {step} or {index}, else one
%   series file. At most two source meshes are held at once. See
%   doc/sequences.md.
%   Options (positional in this order, or name/value):
%     'method'        'linear' (blend the bracketing steps), 'nearest' or
%                     'previous' ('linear')
%     'clamp'         a time outside the source range takes the end step
%                     instead of failing (false)
%     'blend_points'  linear resampling also blends point coordinates
%                     (false)
%     'out_format'    forced output format, [] to infer

    o = mio.internal.opts(varargin, {'method', 'linear', 'clamp', false, ...
        'blend_points', false, 'out_format', []}, 'mio.sequence_resample');
    m = [];
    if ischar(o.method)
        m = find(strcmp(o.method, {'linear', 'nearest', 'previous'}), 1);
    end
    if isempty(m)
        if ischar(o.method)
            name = o.method;
        else
            name = '?';
        end
        error('meshioplusplus:argument', 'meshio++: unknown resample method ''%s''', name);
    end
    meshioplusplus_mex('mio_sequence_resample', seq, mio.internal.str(out_path, 'out_path'), ...
        mio.internal.str(o.out_format, 'out_format'), double(times), m - 1, ...
        logical(o.clamp), logical(o.blend_points));
end
