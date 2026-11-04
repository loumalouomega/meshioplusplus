function exodus_series_write_data(series, time, mesh)
%EXODUS_SERIES_WRITE_DATA  Append one time step to an Exodus series.
%   mio.exodus_series_write_data(SERIES, TIME, MESH) appends a step at the
%   finite time value TIME holding MESH's point and cell data only, so a
%   solver can pass the very mesh it updates in place. Its cell blocks must
%   match those of the grid written by mio.exodus_series_write_points_cells.
%   The first step fixes the field names, dtypes and shapes.

    meshioplusplus_mex('mio_exodus_series_write_data', series, double(time), mesh);
end
