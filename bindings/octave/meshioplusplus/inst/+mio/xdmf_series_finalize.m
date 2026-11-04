function xdmf_series_finalize(series)
%XDMF_SERIES_FINALIZE  Finish an XDMF series. Idempotent.
%   mio.xdmf_series_finalize(SERIES) writes the closing data and closes the
%   file. Releasing the series (or dropping its last reference) also does
%   this, but cannot report a write failure; call this to see one.

    meshioplusplus_mex('mio_xdmf_series_finalize', series);
end
