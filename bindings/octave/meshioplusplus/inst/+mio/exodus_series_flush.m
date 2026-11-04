function exodus_series_flush(series)
%EXODUS_SERIES_FLUSH  Publish the steps written so far.
%   mio.exodus_series_flush(SERIES) makes the file readable up to the last
%   completed step while the series stays open.

    meshioplusplus_mex('mio_exodus_series_flush', series);
end
