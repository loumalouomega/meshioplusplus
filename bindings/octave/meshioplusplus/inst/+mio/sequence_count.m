function n = sequence_count(seq)
%SEQUENCE_COUNT  Number of steps in a sequence.
%   N = mio.sequence_count(SEQ) for a mio.SequenceHandle.

    n = meshioplusplus_mex('mio_sequence_count', seq);
end
