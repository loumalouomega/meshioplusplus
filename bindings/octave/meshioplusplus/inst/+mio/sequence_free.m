function sequence_free(seq)
%SEQUENCE_FREE  Release a sequence handle now.
%   mio.sequence_free(SEQ) frees the native sequence behind SEQ. Idempotent,
%   and a no-op for anything that is not a mio.SequenceHandle.

    meshioplusplus_mex('mio_sequence_free', seq);
end
