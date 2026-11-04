function k = sequence_step(seq, index)
%SEQUENCE_STEP  Step index of entry INDEX within its own file.
%   K = mio.sequence_step(SEQ, INDEX) is the 0-based step of entry INDEX
%   (1-based) inside its file -- the time_step of mio.read -- and 0 for a
%   single-step file.

    k = meshioplusplus_mex('mio_sequence_step', seq, double(index));
end
