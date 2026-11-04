function tf = exodus_series_finalized(series)
%EXODUS_SERIES_FINALIZED  True once an Exodus series has been finalized.

    tf = meshioplusplus_mex('mio_exodus_series_finalized', series);
end
