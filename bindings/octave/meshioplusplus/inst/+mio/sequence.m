function seq = sequence(pattern, varargin)
%SEQUENCE  Open a multi-file / transient sequence from a glob pattern.
%   SEQ = mio.sequence(PATTERN, ...) opens an ordered plan over the files
%   matching PATTERN (or the steps inside one multi-step file), read one
%   step at a time. The handle stores paths, per-file step indices and time
%   values, never a mesh -- which is why mio.sequence_read hands back an
%   independent mesh rather than a borrow the sequence would have to cache.
%   That is what keeps a 500-step dataset traversable.
%
%   Ordering is natural-numeric, so out_9.vtu precedes out_10.vtu.
%   PATTERN is a glob: * and ? only (no **, no [set]); the directory part is
%   taken literally. A pattern matching nothing is an error.
%   Options (positional in this order, or name/value):
%     'format'     forced input format, [] to resolve per file
%     'times'      explicit per-entry times ([]); the count must match
%     'time_from'  'auto', 'file', 'filename' or 'index' ('auto')
%   SEQ is a mio.SequenceHandle; free it with mio.sequence_free or by
%   clearing it.

    o = mio.internal.opts(varargin, {'format', [], 'times', [], 'time_from', 'auto'}, ...
        'mio.sequence');
    seq = meshioplusplus_mex('mio_sequence_open', mio.internal.str(pattern, 'pattern'), ...
        mio.internal.str(o.format, 'format'), mio.internal.seq_times(o.times), ...
        mio.internal.str(o.time_from, 'time_from'), false);
end
