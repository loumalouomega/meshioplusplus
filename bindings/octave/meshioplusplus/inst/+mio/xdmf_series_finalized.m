function tf = xdmf_series_finalized(series)
%XDMF_SERIES_FINALIZED  True once an XDMF series has been finalized.

    tf = meshioplusplus_mex('mio_xdmf_series_finalized', series);
end
