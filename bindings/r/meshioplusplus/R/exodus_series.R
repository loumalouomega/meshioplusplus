#' Stateful Exodus output
#'
#' Fixed geometry, Point/Side sets and `exodus:attr:*` fields are written once.
#' The first data step fixes field names, dtypes and shapes. A netCDF-enabled
#' meshio++ library is required. Flush publishes completed steps; finalize
#' closes the file. Release is idempotent and GC also releases the handle.
#' @param path Output Exodus path.
#' @param series An owning `mio_exodus_series` handle.
#' @param mesh A `mio_mesh` with the fixed grid or a step's fields.
#' @param time Finite simulation time.
#' @return The constructor returns an owning series. Observers return step
#'   count or logical state; mutations return invisibly.
#' @export
mio_exodus_series <- function(path) .Call(R_mio_exodus_series_create, as.character(path))

#' @rdname mio_exodus_series
#' @export
mio_exodus_series_write_points_cells <- function(series, mesh) invisible(.Call(R_mio_exodus_series_write_points_cells, series, mesh))
#' @rdname mio_exodus_series
#' @export
mio_exodus_series_write_data <- function(series, time, mesh) invisible(.Call(R_mio_exodus_series_write_data, series, as.numeric(time), mesh))
#' @rdname mio_exodus_series
#' @export
mio_exodus_series_flush <- function(series) invisible(.Call(R_mio_exodus_series_flush, series))
#' @rdname mio_exodus_series
#' @export
mio_exodus_series_finalize <- function(series) invisible(.Call(R_mio_exodus_series_finalize, series))
#' @rdname mio_exodus_series
#' @export
mio_exodus_series_num_steps <- function(series) .Call(R_mio_exodus_series_num_steps, series)
#' @rdname mio_exodus_series
#' @export
mio_exodus_series_finalized <- function(series) .Call(R_mio_exodus_series_finalized, series)
#' @rdname mio_exodus_series
#' @export
mio_exodus_series_release <- function(series) invisible(.Call(R_mio_exodus_series_release, series))
#' @rdname mio_exodus_series
#' @export
mio_exodus_series_is_open <- function(series) .Call(R_mio_exodus_series_is_open, series)
