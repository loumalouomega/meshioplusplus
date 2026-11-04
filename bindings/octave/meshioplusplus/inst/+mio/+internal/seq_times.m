function t = seq_times(times)
%SEQ_TIMES  The 'times' option of mio.sequence / mio.sequence_list: []
%   (derive the times) or a numeric vector, as double.

    if isempty(times) && ~iscell(times)
        t = [];
        return;
    end
    if ~(isnumeric(times) || islogical(times))
        error('meshioplusplus:argument', '`times` must be a numeric vector or []');
    end
    t = double(times(:));
end
