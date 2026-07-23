# Error types. No `mio_status` code is ever returned to the caller: every
# wrapper that sees a non-zero status or a NULL return throws, carrying the
# message from the C API's thread-local `mio_last_error()`.

"""
    MeshioError(status, msg)

Raised when a meshio++ C API call fails. `status` is the raw `mio_status`
code (see [`STATUS_NAMES`](@ref)) and `msg` is the message the C API recorded
in `mio_last_error()` for the failing call.
"""
struct MeshioError <: Exception
    status::Int
    msg::String
end

"""
Human-readable names of the `mio_status` codes, indexed by `status + 1`.
"""
const STATUS_NAMES = ("ok", "read error", "write error", "invalid argument",
                      "not found", "unsupported", "internal error")

function Base.showerror(io::IO, e::MeshioError)
    name = 0 <= e.status < length(STATUS_NAMES) ? STATUS_NAMES[e.status+1] : "status $(e.status)"
    print(io, "MeshioError (", name, "): ", e.msg)
end

"""
    BorrowError(msg)

Raised when a zero-copy borrow (see [`points_ptr`](@ref),
[`connectivity_ptr`](@ref)) is used outside its validity window — that is,
after a mutating call on the owning [`Mesh`](@ref) or after the mesh was
closed. The C API's rule is that a borrowed pointer stays valid only until the
next mutating `mio_mesh_*` call on that mesh; this exception is what turns a
violation into a clean error instead of a read of freed or reallocated memory.
"""
struct BorrowError <: Exception
    msg::String
end

Base.showerror(io::IO, e::BorrowError) = print(io, "BorrowError: ", e.msg)

# Throwing helpers -----------------------------------------------------------

"""Throw a `MeshioError` for `status`, reading the C API's last message."""
@noinline function _throw_status(status::Integer)
    throw(MeshioError(Int(status), last_error()))
end

"""Check a `mio_status` return value; throw unless it is `MIO_OK` (0)."""
function _check(status::Integer)
    status == 0 || _throw_status(status)
    nothing
end

"""Check a handle-returning call; throw when it returned NULL."""
function _check_ptr(ptr::Ptr{Cvoid})
    ptr == C_NULL && throw(MeshioError(MIO_ERR_INTERNAL, last_error()))
    ptr
end

"""Check a `-1`-on-error integer return value."""
function _check_count(n::Integer, what::AbstractString)
    n < 0 && throw(MeshioError(MIO_ERR_INVALID_ARG, isempty(last_error()) ?
                               "$what failed" : last_error()))
    n
end
