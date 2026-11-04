function t = sequence_time(seq, index)
%SEQUENCE_TIME  Time value of entry INDEX (1-based) of a sequence.
%   T = mio.sequence_time(SEQ, INDEX).

    t = meshioplusplus_mex('mio_sequence_time', seq, double(index));
end
