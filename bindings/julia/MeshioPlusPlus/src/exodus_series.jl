"""
    ExodusSeries(path)
    ExodusSeries(f::Function, path)

Stateful Exodus writer (requires netCDF). `write_points_cells!` fixes geometry,
sets and `exodus:attr:*` fields. `write_data!` appends fields with a stable schema.
`flush!` publishes completed steps; `finalize!` closes the file. `close` releases
the owning handle and is idempotent. The do-block form always closes the handle.
"""
mutable struct ExodusSeries
    ptr::Ptr{Cvoid}
    function ExodusSeries(ptr::Ptr{Cvoid})
        s = new(ptr)
        finalizer(s) do x
            if x.ptr != C_NULL
                ccall(_sym(:mio_exodus_series_free), Cvoid, (Ptr{Cvoid},), x.ptr)
                x.ptr = C_NULL
            end
        end
        s
    end
end

ExodusSeries(path::AbstractString) = ExodusSeries(_check_ptr(ccall(
    _sym(:mio_exodus_series_create), Ptr{Cvoid}, (Cstring,), path)))

function ExodusSeries(f::Function, path::AbstractString)
    s = ExodusSeries(path)
    try
        f(s)
    finally
        close(s)
    end
end

function _handle(s::ExodusSeries)
    s.ptr == C_NULL && throw(MeshioError(MIO_ERR_INVALID_ARG, "Exodus series has been closed"))
    s.ptr
end

function Base.close(s::ExodusSeries)
    if s.ptr != C_NULL
        ccall(_sym(:mio_exodus_series_free), Cvoid, (Ptr{Cvoid},), s.ptr)
        s.ptr = C_NULL
    end
    nothing
end

Base.isopen(s::ExodusSeries) = s.ptr != C_NULL
write_points_cells!(s::ExodusSeries, m::Mesh) = _check(ccall(
    _sym(:mio_exodus_series_write_points_cells), Cint,
    (Ptr{Cvoid}, Ptr{Cvoid}), _handle(s), _handle(m)))
write_data!(s::ExodusSeries, time::Real, m::Mesh) = _check(ccall(
    _sym(:mio_exodus_series_write_data), Cint,
    (Ptr{Cvoid}, Cdouble, Ptr{Cvoid}), _handle(s), Cdouble(time), _handle(m)))
flush!(s::ExodusSeries) = _check(ccall(
    _sym(:mio_exodus_series_flush), Cint, (Ptr{Cvoid},), _handle(s)))
finalize!(s::ExodusSeries) = _check(ccall(
    _sym(:mio_exodus_series_finalize), Cint, (Ptr{Cvoid},), _handle(s)))
finalized(s::ExodusSeries) = ccall(
    _sym(:mio_exodus_series_finalized), Int32, (Ptr{Cvoid},), _handle(s)) == 1
function num_steps(s::ExodusSeries)
    n = ccall(_sym(:mio_exodus_series_num_steps), Int64, (Ptr{Cvoid},), _handle(s))
    n < 0 && throw(MeshioError(MIO_ERR_INVALID_ARG, "cannot query Exodus series steps"))
    n
end
