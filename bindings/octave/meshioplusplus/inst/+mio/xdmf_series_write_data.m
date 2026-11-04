function xdmf_series_write_data(series, time, mesh)
%XDMF_SERIES_WRITE_DATA  Append one time step to an XDMF series.
%   mio.xdmf_series_write_data(SERIES, TIME, MESH) appends a step at the
%   time value TIME holding MESH's point and cell data only, so a
%   solver can pass the very mesh it updates in place. Its cell blocks must
%   match those of the grid written by mio.xdmf_series_write_points_cells.

    meshioplusplus_mex('mio_xdmf_series_write_data', series, double(time), mesh);
end
