function exodus_series_finalize(series)
%EXODUS_SERIES_FINALIZE  Finish an Exodus series. Idempotent.
%   mio.exodus_series_finalize(SERIES) writes the closing data and closes the
%   file. Releasing the series (or dropping its last reference) also does
%   this, but cannot report a write failure; call this to see one.

    meshioplusplus_mex('mio_exodus_series_finalize', series);
end
