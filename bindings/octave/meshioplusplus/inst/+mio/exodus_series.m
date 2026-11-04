function s = exodus_series(path)
%EXODUS_SERIES  Open a stateful Exodus II series writer.
%   S = mio.exodus_series(PATH) starts an Exodus file at PATH and returns a
%   mio.ExodusSeries. Fixed geometry, Point/Side sets and exodus:attr:*
%   fields are written once with mio.exodus_series_write_points_cells; each
%   step's fields are appended with mio.exodus_series_write_data, and the
%   first step fixes field names, dtypes and shapes.
%
%   A netCDF-enabled meshio++ library is required: without it construction
%   fails with an error naming MESHIOPLUSPLUS_WITH_NETCDF.
%   mio.exodus_series_flush publishes completed steps;
%   mio.exodus_series_finalize closes the file. mio.exodus_series_release is
%   idempotent, and dropping the last reference also releases the handle.
%
%   Lifecycle: mio.exodus_series_num_steps, mio.exodus_series_finalized,
%   mio.exodus_series_is_open.

    s = meshioplusplus_mex('mio_exodus_series_create', mio.internal.str(path, 'path'));
end
