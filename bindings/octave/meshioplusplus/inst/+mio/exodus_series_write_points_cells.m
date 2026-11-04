function exodus_series_write_points_cells(series, mesh)
%EXODUS_SERIES_WRITE_POINTS_CELLS  Write the fixed grid of an Exodus series.
%   mio.exodus_series_write_points_cells(SERIES, MESH) writes MESH's fixed
%   geometry, Point/Side sets and exodus:attr:* fields once, before any
%   step. SERIES is a mio.ExodusSeries.

    meshioplusplus_mex('mio_exodus_series_write_points_cells', series, mesh);
end
