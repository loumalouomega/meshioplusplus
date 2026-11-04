function n = exodus_series_num_steps(series)
%EXODUS_SERIES_NUM_STEPS  Number of time steps written so far (a double).

    n = meshioplusplus_mex('mio_exodus_series_num_steps', series);
end
