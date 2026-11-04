function xdmf_series_flush(series)
%XDMF_SERIES_FLUSH  Publish the steps written so far.
%   mio.xdmf_series_flush(SERIES) makes the file readable up to the last
%   completed step while the series stays open.
%   It re-serializes the whole .xdmf, so flush on your own cadence rather
%   than after every step.

    meshioplusplus_mex('mio_xdmf_series_flush', series);
end
