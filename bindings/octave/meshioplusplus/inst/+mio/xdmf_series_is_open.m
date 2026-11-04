function tf = xdmf_series_is_open(series)
%XDMF_SERIES_IS_OPEN  True while SERIES is a live mio.XdmfSeries.
%   False once released, and for anything that is not a mio.XdmfSeries.

    tf = meshioplusplus_mex('mio_xdmf_series_is_open', series);
end
