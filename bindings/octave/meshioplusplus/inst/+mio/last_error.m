function s = last_error()
%LAST_ERROR  The message of the most recent failed call on this thread, or ''.

    s = meshioplusplus_mex('mio_last_error');
end
