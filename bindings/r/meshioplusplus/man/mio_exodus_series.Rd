\name{mio_exodus_series}
\alias{mio_exodus_series}
\alias{mio_exodus_series_write_points_cells}
\alias{mio_exodus_series_write_data}
\alias{mio_exodus_series_flush}
\alias{mio_exodus_series_finalize}
\alias{mio_exodus_series_num_steps}
\alias{mio_exodus_series_finalized}
\alias{mio_exodus_series_release}
\alias{mio_exodus_series_is_open}
\title{Stateful Exodus output}
\usage{
mio_exodus_series(path)
mio_exodus_series_write_points_cells(series, mesh)
mio_exodus_series_write_data(series, time, mesh)
mio_exodus_series_flush(series)
mio_exodus_series_finalize(series)
mio_exodus_series_num_steps(series)
mio_exodus_series_finalized(series)
mio_exodus_series_release(series)
mio_exodus_series_is_open(series)
}
\arguments{
\item{path}{Output Exodus path.}
\item{series}{An owning \code{mio_exodus_series} handle.}
\item{mesh}{A \code{mio_mesh} containing the fixed grid or a step's fields.}
\item{time}{Finite simulation time.}
}
\value{The constructor returns an owning series. Observers return step count or logical state; mutations return invisibly.}
\description{Requires a netCDF-enabled library. Write fixed geometry, Point/Side sets and \code{exodus:attr:*} arrays once, then append steps with stable field names, dtypes and shapes. Flush publishes completed steps. Finalize closes the file; release is idempotent and GC also releases the handle.}
