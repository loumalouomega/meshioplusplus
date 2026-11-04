function s = sequence_time_source(seq, index)
%SEQUENCE_TIME_SOURCE  Where entry INDEX's time came from.
%   S = mio.sequence_time_source(SEQ, INDEX) is 'explicit', 'file',
%   'filename' or 'index' (the fallback). Reported rather than left to be
%   guessed: "the file said 0.25" and "nothing said anything, so this is
%   position 3" are different facts.

    s = meshioplusplus_mex('mio_sequence_time_source', seq, double(index));
end
