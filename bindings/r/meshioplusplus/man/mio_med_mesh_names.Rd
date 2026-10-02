\name{mio_med_mesh_names}
\alias{mio_med_mesh_names}
\alias{mio_med_read_named}
\alias{mio_med_write_multi}
\title{Named MED meshes}
\usage{
mio_med_mesh_names(path)
mio_med_read_named(path, name, time_step = 0L, lenient = FALSE)
mio_med_write_multi(path, meshes, names, version = "4.1.0")
}
\arguments{
\item{path}{MED file path.}
\item{name}{A mesh name.}
\item{time_step}{Zero-based step index; negative counts from the end.}
\item{lenient}{Skip unsupported field constructs.}
\item{meshes}{Nonempty list of mesh handles.}
\item{names}{Unique, nonempty names, one per mesh.}
\item{version}{MED version string, default \code{"4.1.0"}.}
}
\value{Names, an owning mesh handle, or invisibly NULL for writes.}
\description{Enumerate or read named MED meshes, or write several meshes to one file. Ordinary named nodal and element profiles expand with NaN fill. Requires an HDF5-enabled native library.}
