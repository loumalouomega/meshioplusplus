# The raw interface to `libmeshioplusplus`: library discovery, symbol lookup,
# the C enums and structs mirrored bit-for-bit, and the small helpers that
# implement the C API's calling conventions (the string-getter protocol, name
# lists, dtype mapping). Everything above this file calls
# `ccall(_sym(:mio_...), ...)` directly rather than through a per-symbol
# forwarder -- the C API is flat, so a forwarder layer would only duplicate it.

# --- library discovery -------------------------------------------------------

const _HANDLE = Ref{Ptr{Cvoid}}(C_NULL)
const _SYMS = Dict{Symbol,Ptr{Cvoid}}()
const _LIBPATH = Ref{String}("")

const _BUILD_HINT = """
meshio++: could not load the C API shared library `libmeshioplusplus`.

This package binds the installed C library; build and install it first:

    ./build/configure.sh --c-api --build
    cmake --install build/cpp-release --prefix /your/prefix

then either put /your/prefix/lib on the loader path, or point this package
straight at the library:

    ENV["MESHIOPLUSPLUS_LIB"] = "/your/prefix/lib/libmeshioplusplus.so"

See the meshio++ docs (doc/julia.md, doc/c_api.md) for details.
"""

"""
    library_path() -> String

Path of the `libmeshioplusplus` shared library this session loaded.
"""
library_path() = _LIBPATH[]

function _load_library()
    explicit = get(ENV, "MESHIOPLUSPLUS_LIB", "")
    if !isempty(explicit)
        handle = try
            Libdl.dlopen(explicit)
        catch err
            error(_BUILD_HINT * "\nMESHIOPLUSPLUS_LIB=\"$explicit\" could not be " *
                  "opened: $(sprint(showerror, err))")
        end
        _HANDLE[] = handle
        _LIBPATH[] = explicit
        return
    end
    found = Libdl.find_library(["libmeshioplusplus", "meshioplusplus"])
    isempty(found) && error(_BUILD_HINT)
    _HANDLE[] = Libdl.dlopen(found)
    _LIBPATH[] = Libdl.dlpath(_HANDLE[])
    return
end

"""Resolve (and memoize) one `mio_*` symbol in the loaded library."""
function _sym(name::Symbol)
    get!(_SYMS, name) do
        _HANDLE[] == C_NULL && error(_BUILD_HINT)
        Libdl.dlsym(_HANDLE[], name)
    end
end

# --- enums (must match bindings/c/include/meshioplusplus/meshioplusplus.h) ---

const MIO_OK = 0
const MIO_ERR_READ = 1
const MIO_ERR_WRITE = 2
const MIO_ERR_INVALID_ARG = 3
const MIO_ERR_NOT_FOUND = 4
const MIO_ERR_UNSUPPORTED = 5
const MIO_ERR_INTERNAL = 6

const MIO_MAX_NDIM = 8

# mio_data_location
const MIO_DATA_POINT = Cint(0)
const MIO_DATA_CELL = Cint(1)
const MIO_DATA_FIELD = Cint(2)

# mio_region_kind
const MIO_REGION_POINT = Cint(0)
const MIO_REGION_CELL = Cint(1)
const MIO_REGION_SIDE = Cint(2)

# mio_diff_verdict
const MIO_DIFF_IDENTICAL = 0
const MIO_DIFF_EQUAL_WITHIN_TOLERANCE = 1
const MIO_DIFF_DIFFERENT = 2

# mio_dtype, in enum order: the Julia element type of each `mio_dtype` value.
const _DTYPES = (Float32, Float64, Int8, Int16, Int32, Int64,
                 UInt8, UInt16, UInt32, UInt64)

"""Julia element type for a `mio_dtype` value."""
function _jltype(dt::Integer)
    0 <= dt < length(_DTYPES) ||
        throw(MeshioError(MIO_ERR_INVALID_ARG, "unknown mio_dtype $dt"))
    _DTYPES[dt+1]
end

"""`mio_dtype` value for a Julia element type."""
function _miodtype(::Type{T}) where {T}
    i = findfirst(==(T), _DTYPES)
    i === nothing && throw(MeshioError(MIO_ERR_INVALID_ARG,
                                       "no meshio++ dtype for Julia type $T"))
    Cint(i - 1)
end

# --- C structs mirrored bit-for-bit -----------------------------------------

"""
Mirror of C `mio_read_opts`. Field order, types and the trailing `reserved`
padding are ABI: they must match `meshioplusplus.h` exactly. Always build one
through [`ReadOptions`](@ref) rather than by hand.
"""
struct _CReadOpts
    points_only::Cint
    metadata_only::Cint
    arrays::Ptr{Cstring}
    num_arrays::Int64
    mmap_mode::Cint
    _pad::Cint
    reserved::NTuple{6,Int64}
end

"""Mirror of C `mio_stats_report`."""
struct _CStatsReport
    num_points::Int64
    num_cells::Int64
    bbox_min::NTuple{3,Cdouble}
    bbox_max::NTuple{3,Cdouble}
    extent::NTuple{3,Cdouble}
    centroid::NTuple{3,Cdouble}
    total_area::Cdouble
    signed_volume::Cdouble
    unsigned_volume::Cdouble
    num_inverted::Int64
end

"""Mirror of C `mio_data_array_info`."""
struct _CDataArrayInfo
    location::Cint
    dtype::Cint
    num_blocks::Int64
    num_entries::Int64
    num_components::Int64
    num_values::Int64
    min::Cdouble
    max::Cdouble
    mean::Cdouble
    num_nan::Int64
    num_inf::Int64
    num_finite::Int64
    inconsistent_blocks::Cint
end

"""Mirror of C `mio_region_info`."""
struct _CRegionInfo
    kind::Cint
    dim::Cint
    tag::Int64
    num_entries::Int64
    stride::Int64
end

# Layout guards: a mismatch here would corrupt every call taking these structs,
# silently. Checked once at load rather than trusted.
function _check_abi_layout()
    sizeof(_CReadOpts) == 80 ||
        error("meshio++: mio_read_opts layout mismatch ($(sizeof(_CReadOpts)) bytes)")
    sizeof(_CStatsReport) == 144 ||
        error("meshio++: mio_stats_report layout mismatch ($(sizeof(_CStatsReport)) bytes)")
    sizeof(_CRegionInfo) == 32 ||
        error("meshio++: mio_region_info layout mismatch ($(sizeof(_CRegionInfo)) bytes)")
    nothing
end

# --- calling-convention helpers ---------------------------------------------

"""
    _getstring(f) -> String

Run the C API's string-getter protocol (rule 5 of the header): call `f(buf,
buflen)` once with a NULL buffer to learn the required length, then again with
a buffer of that size. `f` must return the length excluding the NUL, or -1.
"""
function _getstring(f)
    n = f(Ptr{UInt8}(C_NULL), Int64(0))
    n < 0 && _throw_status(MIO_ERR_INVALID_ARG)
    n == 0 && return ""
    buf = Vector{UInt8}(undef, n + 1)
    m = GC.@preserve buf f(pointer(buf), Int64(n + 1))
    m < 0 && _throw_status(MIO_ERR_INVALID_ARG)
    String(buf[1:min(n, m)])
end

"""
    _with_names(f, names) -> Any

Marshal a `Vector{<:AbstractString}` into the C API's `(const char* const*,
int64_t count)` name-list convention and call `f(ptr, count)`. The Julia
strings stay rooted for the duration of the call.
"""
function _with_names(f, names)
    isempty(names) && return f(Ptr{Cstring}(C_NULL), Int64(0))
    roots = [Base.cconvert(Cstring, String(n)) for n in names]
    ptrs = Cstring[Base.unsafe_convert(Cstring, r) for r in roots]
    GC.@preserve roots ptrs begin
        f(pointer(ptrs), Int64(length(ptrs)))
    end
end
