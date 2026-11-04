function exodus_series_release(series)
%EXODUS_SERIES_RELEASE  Free an Exodus series now. Idempotent.
%   mio.exodus_series_release(SERIES) finalizes the file if it was not finalized and
%   frees the native writer; a second call is a no-op, and any other call on
%   SERIES then raises. A non-series argument is an error.

    meshioplusplus_mex('mio_exodus_series_release', series);
end
