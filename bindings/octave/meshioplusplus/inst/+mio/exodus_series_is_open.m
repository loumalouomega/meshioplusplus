function tf = exodus_series_is_open(series)
%EXODUS_SERIES_IS_OPEN  True while SERIES is a live mio.ExodusSeries.
%   False once released, and for anything that is not a mio.ExodusSeries.

    tf = meshioplusplus_mex('mio_exodus_series_is_open', series);
end
