function seq = sequence_list(paths, varargin)
%SEQUENCE_LIST  Open a sequence from an explicit, ordered path list.
%   SEQ = mio.sequence_list(PATHS, ...) opens a sequence over the cellstr
%   PATHS. The order is yours and is kept unless 'sort' is true.
%   Options (positional in this order, or name/value):
%     'format'     forced input format, [] to resolve per file
%     'times'      explicit per-entry times ([])
%     'time_from'  'auto', 'file', 'filename' or 'index' ('auto')
%     'sort'       natural-numeric sort the list too (false)
%   SEQ is a mio.SequenceHandle.

    o = mio.internal.opts(varargin, {'format', [], 'times', [], 'time_from', 'auto', ...
        'sort', false}, 'mio.sequence_list');
    seq = meshioplusplus_mex('mio_sequence_open_list', mio.internal.cellstr_arg(paths, 'paths'), ...
        mio.internal.str(o.format, 'format'), mio.internal.seq_times(o.times), ...
        mio.internal.str(o.time_from, 'time_from'), logical(o.sort));
end
