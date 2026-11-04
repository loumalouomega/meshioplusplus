function xdmf_series_write_points_cells(series, mesh)
%XDMF_SERIES_WRITE_POINTS_CELLS  Write the fixed grid of an XDMF series.
%   mio.xdmf_series_write_points_cells(SERIES, MESH) writes MESH's points,
%   cells and fixed regions once, before any step. SERIES is a
%   mio.XdmfSeries.

    meshioplusplus_mex('mio_xdmf_series_write_points_cells', series, mesh);
end
