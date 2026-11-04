function m = sequence_read(seq, index)
%SEQUENCE_READ  Read one step of a sequence.
%   M = mio.sequence_read(SEQ, INDEX) reads entry INDEX (1-based) as a
%   mio.MeshHandle, independent of the sequence, which caches nothing.

    m = meshioplusplus_mex('mio_sequence_read', seq, double(index));
end
