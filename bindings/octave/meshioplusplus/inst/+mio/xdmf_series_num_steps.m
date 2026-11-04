function n = xdmf_series_num_steps(series)
%XDMF_SERIES_NUM_STEPS  Number of time steps written so far (a double).

    n = meshioplusplus_mex('mio_xdmf_series_num_steps', series);
end
