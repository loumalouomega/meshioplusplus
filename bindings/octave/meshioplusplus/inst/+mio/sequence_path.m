function p = sequence_path(seq, index)
%SEQUENCE_PATH  File path of entry INDEX (1-based) of a sequence.
%   P = mio.sequence_path(SEQ, INDEX).

    p = meshioplusplus_mex('mio_sequence_path', seq, double(index));
end
