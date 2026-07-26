# The transient (time-series) XDMF writer: the write half of what
# `ReadOptions(time_step=...)` and `read_metadata(...).time_values` expose on
# the read side.
#
# This is the one writer in meshio++ that `write` cannot express. A series is a
# stateful multi-call object -- the grid goes out once and each solve appends a
# cheap step -- so it has no `(path, mesh)` shape for the format registry to
# name, and it gets its own handle here rather than a keyword on `write`.

"""
    XdmfSeries(path; data_format="HDF", gzip_level=-1)
    XdmfSeries(f::Function, path; data_format="HDF", gzip_level=-1)

A transient XDMF writer. The mesh is written **once** with
[`write_points_cells!`](@ref) and each solve appends one step with
[`write_data!`](@ref):

```julia
s = XdmfSeries("simulation.xdmf")
write_points_cells!(s, mesh)
for k in 0:nsteps-1
    solve!(mesh)
    write_data!(s, k * dt, mesh)
end
close(s)
```

`data_format` is `"HDF"` (the default; needs an HDF5-enabled library), `"XML"`
(everything inline in the `.xdmf`) or `"Binary"`. `gzip_level` applies to
`"HDF"` datasets only; negative means no compression. An unknown format, or
`"HDF"` against a library built without HDF5, throws a [`MeshioError`](@ref)
from the constructor.

The `.xdmf` light data is **buffered until the series is finalized**, so the
file is only readable after [`finalize!`](@ref) or [`close`](@ref). Heavy data
for `"HDF"` goes to a `<path minus extension>.h5` *sibling* of the `.xdmf`.

Like [`Mesh`](@ref), the handle is released by a finalizer; [`close`](@ref)
does it deterministically and is idempotent. The function form runs `f(series)`
and closes the series afterwards even if `f` throws:

```julia
XdmfSeries("simulation.xdmf") do s
    write_points_cells!(s, mesh)
    write_data!(s, 0.0, mesh)
end
```
"""
mutable struct XdmfSeries
    ptr::Ptr{Cvoid}

    function XdmfSeries(ptr::Ptr{Cvoid})
        s = new(ptr)
        finalizer(s) do x
            if x.ptr != C_NULL
                ccall(_sym(:mio_xdmf_series_free), Cvoid, (Ptr{Cvoid},), x.ptr)
                x.ptr = C_NULL
            end
        end
        s
    end
end

function XdmfSeries(path::AbstractString; data_format::AbstractString="HDF",
                    gzip_level::Integer=-1)
    XdmfSeries(_check_ptr(ccall(_sym(:mio_xdmf_series_create), Ptr{Cvoid},
                                (Cstring, Cstring, Int32),
                                path, data_format, Int32(gzip_level))))
end

function XdmfSeries(f::Function, path::AbstractString; kwargs...)
    s = XdmfSeries(path; kwargs...)
    try
        f(s)
    finally
        close(s)
    end
end

"""Raw handle of an open series; throws once it has been closed."""
function _handle(s::XdmfSeries)
    s.ptr == C_NULL &&
        throw(MeshioError(MIO_ERR_INVALID_ARG, "xdmf series has been closed"))
    s.ptr
end

"""
    close(series::XdmfSeries)

Finalize (writing the `.xdmf`) and release the handle immediately instead of
waiting for the finalizer. Idempotent.

A write failure during this implicit finalize is swallowed, exactly as in the C
API — call [`finalize!`](@ref) first to see one.
"""
function Base.close(s::XdmfSeries)
    if s.ptr != C_NULL
        ccall(_sym(:mio_xdmf_series_free), Cvoid, (Ptr{Cvoid},), s.ptr)
        s.ptr = C_NULL
    end
    nothing
end

"""    isopen(series::XdmfSeries) -> Bool"""
Base.isopen(s::XdmfSeries) = s.ptr != C_NULL

function Base.show(io::IO, s::XdmfSeries)
    if !isopen(s)
        print(io, "XdmfSeries(<closed>)")
    else
        print(io, "XdmfSeries(", num_steps(s), " steps written)")
    end
end

"""
    write_points_cells!(series, mesh)

Write the static grid every step shares. Call once, before the first
[`write_data!`](@ref). Only the mesh's points and cells are used.
"""
function write_points_cells!(s::XdmfSeries, m::Mesh)
    _check(ccall(_sym(:mio_xdmf_series_write_points_cells), Cint,
                 (Ptr{Cvoid}, Ptr{Cvoid}), _handle(s), _handle(m)))
end

"""
    write_data!(series, time, mesh)

Write one step's `point_data` and `cell_data` at `time`. The mesh's geometry is
ignored — a solver can pass the very object it is updating in place — but its
cell blocks must match those of the static grid.
"""
function write_data!(s::XdmfSeries, time::Real, m::Mesh)
    _check(ccall(_sym(:mio_xdmf_series_write_data), Cint,
                 (Ptr{Cvoid}, Cdouble, Ptr{Cvoid}),
                 _handle(s), Cdouble(time), _handle(m)))
end

"""
    finalize!(series)

Write the `.xdmf` and close the heavy-data container. Idempotent, and
[`close`](@ref) does it too — call this explicitly so a write failure surfaces
as a [`MeshioError`](@ref) instead of being swallowed.

Named `finalize!` rather than `finalize` to stay clear of `Base.finalize`,
which runs an object's GC finalizer and means something quite different.
"""
function finalize!(s::XdmfSeries)
    _check(ccall(_sym(:mio_xdmf_series_finalize), Cint, (Ptr{Cvoid},), _handle(s)))
end

"""
    num_steps(series) -> Int

Number of steps written so far.
"""
function num_steps(s::XdmfSeries)
    Int(_check_count(ccall(_sym(:mio_xdmf_series_num_steps), Int64, (Ptr{Cvoid},),
                           _handle(s)), "num_steps"))
end
