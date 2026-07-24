# The `Mesh` handle, the borrow window, file I/O, and the setter/getter pairs.
#
# Two conventions run through this whole file (see doc/julia.md):
#
#   * COLUMN-MAJOR. Julia is column-major and the C core row-major, so points
#     shaped `(dim, num_points)` and connectivity `(nodes_per_cell, num_cells)`
#     are the SAME MEMORY as the C API's `(num_points, dim)` /
#     `(num_cells, nodes_per_cell)`. Nothing is ever transposed.
#   * 1-BASED. Connectivity is 0-based at the ABI. The +-1 shift happens inside
#     the COPYING accessors only; the zero-copy borrows are left un-shifted,
#     because shifting a borrow would mean copying it.

# --- the handle --------------------------------------------------------------

"""
    Mesh()

An empty meshio++ mesh. Wraps the C API's opaque `mio_mesh*`.

The handle is released by a **finalizer** calling `mio_mesh_free` — unlike the
Fortran module, where meshes are freed explicitly with `call m%free()`. Call
[`close`](@ref) to release one deterministically; it is idempotent and the
finalizer then does nothing.

Every `Mesh` **owns** its handle. Operations that return an opaque C result
(split, partition, reorder, refine, decimate, convert_cells) always transfer
ownership out of the result rather than handing back a borrow into it, so a
mesh can never dangle when its producing result is freed.
"""
mutable struct Mesh
    ptr::Ptr{Cvoid}
    gen::Int

    function Mesh(ptr::Ptr{Cvoid})
        m = new(ptr, 0)
        finalizer(m) do x
            if x.ptr != C_NULL
                ccall(_sym(:mio_mesh_free), Cvoid, (Ptr{Cvoid},), x.ptr)
                x.ptr = C_NULL
            end
        end
        m
    end
end

Mesh() = Mesh(_check_ptr(ccall(_sym(:mio_mesh_create), Ptr{Cvoid}, ())))

"""Raw handle of an open mesh; throws when the mesh has been closed."""
function _handle(m::Mesh)
    m.ptr == C_NULL && throw(MeshioError(MIO_ERR_INVALID_ARG, "mesh has been closed"))
    m.ptr
end

"""
    close(mesh::Mesh)

Release the mesh's C handle immediately instead of waiting for the finalizer.
Idempotent. Every outstanding borrow into the mesh is invalidated, and using
one afterwards raises [`BorrowError`](@ref) rather than reading freed memory.
"""
function Base.close(m::Mesh)
    if m.ptr != C_NULL
        ccall(_sym(:mio_mesh_free), Cvoid, (Ptr{Cvoid},), m.ptr)
        m.ptr = C_NULL
    end
    m.gen += 1  # kill every outstanding borrow
    nothing
end

"""    isopen(mesh::Mesh) -> Bool"""
Base.isopen(m::Mesh) = m.ptr != C_NULL

function Base.show(io::IO, m::Mesh)
    if !isopen(m)
        print(io, "Mesh(<closed>)")
    else
        print(io, "Mesh(", num_points(m), " points, ", num_cells(m), " cells in ",
              num_cell_blocks(m), " blocks)")
    end
end

"""Mark the mesh mutated: every borrow taken before this call becomes invalid."""
_touch!(m::Mesh) = (m.gen += 1; nothing)

# --- the borrow window -------------------------------------------------------

"""
    MeshBorrow{T,N} <: AbstractArray{T,N}

A zero-copy view of mesh-owned memory, obtained from one of the `*_ptr`
accessors. No data is copied: the array aliases the buffer inside the C++ core.

Two things make that safe. The borrow holds a reference to its owning
[`Mesh`](@ref), so the mesh cannot be garbage-collected while the view is
alive; and it records the mesh's mutation generation, so using it after a
mutating call (`set_points!`, `add_cell_block!`, `add_point_data!`,
`append_cell_data!`, `add_field_data!`, `add_region!`) or after
[`close`](@ref) raises [`BorrowError`](@ref).

This mirrors the C API's rule 3: *returned data pointers alias mesh-owned
memory and stay valid until the next mutating `mio_mesh_*` call on that mesh
or `mio_mesh_free()` — read-only accessors never invalidate them.*

`parent(b)` returns the underlying `Array` for hot loops. That escape hatch is
**unchecked**: it is valid only for as long as the borrow itself is.
"""
struct MeshBorrow{T,N} <: AbstractArray{T,N}
    data::Array{T,N}
    owner::Mesh
    gen::Int
end

@inline function _check_borrow(b::MeshBorrow)
    if b.owner.ptr == C_NULL
        throw(BorrowError("the mesh this borrow views has been closed"))
    elseif b.owner.gen != b.gen
        throw(BorrowError("the mesh was mutated after this borrow was taken; " *
                          "take a fresh one (or use a copying accessor)"))
    end
    nothing
end

Base.size(b::MeshBorrow) = (_check_borrow(b); size(b.data))
Base.IndexStyle(::Type{<:MeshBorrow}) = IndexLinear()
Base.@propagate_inbounds function Base.getindex(b::MeshBorrow, i::Int)
    _check_borrow(b)
    b.data[i]
end
Base.@propagate_inbounds function Base.setindex!(b::MeshBorrow, v, i::Int)
    _check_borrow(b)
    b.data[i] = v
end

"""
    parent(b::MeshBorrow) -> Array

The raw aliasing array, without the validity check. Unsafe outside the borrow
window; prefer indexing the `MeshBorrow` itself unless the check shows up in a
profile.
"""
Base.parent(b::MeshBorrow) = b.data

"""Wrap `ptr` as a borrow of `dims` elements of type `T` owned by `m`."""
function _borrow(m::Mesh, ptr::Ptr{Cvoid}, ::Type{T}, dims::Tuple) where {T}
    data = unsafe_wrap(Array, Ptr{T}(ptr), dims; own=false)
    MeshBorrow{T,length(dims)}(data, m, m.gen)
end

# --- version / build introspection ------------------------------------------

"""    version() -> String — the meshio++ version of the loaded C library."""
version() = unsafe_string(ccall(_sym(:mio_version), Cstring, ()))

"""    mesh_backend() -> String — the compile-time mesh backend ("meshio", "native", "kratos")."""
mesh_backend() = unsafe_string(ccall(_sym(:mio_mesh_backend), Cstring, ()))

"""
    last_error() -> String

The message of the most recent failed C API call on this thread ("" if none).
Wrappers surface this through [`MeshioError`](@ref); it is exported for
diagnostics.
"""
last_error() = unsafe_string(ccall(_sym(:mio_last_error), Cstring, ()))

"""    format_readable(format) -> Bool — is `format` readable in this build?"""
format_readable(format::AbstractString) =
    ccall(_sym(:mio_format_readable), Cint, (Cstring,), format) != 0

"""    format_writable(format) -> Bool — is `format` writable in this build?"""
format_writable(format::AbstractString) =
    ccall(_sym(:mio_format_writable), Cint, (Cstring,), format) != 0

"""
    sniff_format(path) -> String

Guess a mesh file's format from its leading bytes. Returns `""` when the
content is not a confident match (the C API never guesses on an ambiguous
magic).
"""
function sniff_format(path::AbstractString)
    _getstring() do buf, len
        ccall(_sym(:mio_sniff_format), Int64, (Cstring, Ptr{UInt8}, Int64), path, buf, len)
    end
end

# --- cell-type metadata ------------------------------------------------------

"""
    cell_type_num_nodes(name) -> Int

Fixed nodes per cell of a meshio cell-type name (e.g. `"tetra10"` → 10), or
`-1` for the variable-size types (polygon, polyhedron, VTK Lagrange) and for
names outside the table.
"""
function cell_type_num_nodes(name::AbstractString)
    t = ccall(_sym(:mio_cell_type_from_name), Cint, (Cstring,), name)
    Int(ccall(_sym(:mio_cell_type_num_nodes), Cint, (Cint,), t))
end

"""
    cell_type_dimension(name) -> Int

Topological dimension (0–3) of a meshio cell-type name, or `-1` when the name
is not in the table.
"""
function cell_type_dimension(name::AbstractString)
    t = ccall(_sym(:mio_cell_type_from_name), Cint, (Cstring,), name)
    Int(ccall(_sym(:mio_cell_type_dimension), Cint, (Cint,), t))
end

# --- file I/O ----------------------------------------------------------------

"""
    ReadOptions(; points_only=false, metadata_only=false, arrays=nothing, mmap=:auto,
                  time_step=0)

Narrow what a read materializes (see `doc/selective_read.md`).

* `points_only` — read geometry only, skipping every data array.
* `metadata_only` — read the header/summary only.
* `arrays` — `nothing` keeps every array; a vector of names keeps only those.
  An **empty vector** means "no arrays at all", which is deliberately distinct
  from `nothing`.
* `mmap` — `:auto`, `:on` or `:off` (see `doc/mmap.md`).

Formats without a native selective path are read in full; the options are
still honoured, just without the saving.
"""
struct ReadOptions
    points_only::Bool
    metadata_only::Bool
    arrays::Union{Nothing,Vector{String}}
    mmap::Symbol

    function ReadOptions(; points_only::Bool=false, metadata_only::Bool=false,
                         arrays=nothing, mmap::Symbol=:auto)
        mmap in (:auto, :on, :off) ||
            throw(ArgumentError("mmap must be :auto, :on or :off, got :$mmap"))
        new(points_only, metadata_only,
            arrays === nothing ? nothing : String[String(a) for a in arrays], mmap)
    end
end

_mmap_mode(s::Symbol) = s === :auto ? Cint(0) : s === :on ? Cint(1) : Cint(2)

"""Build the C struct for `opts` and run `f(Ref{_CReadOpts})` with the names rooted."""
function _with_read_opts(f, opts::ReadOptions)
    base = Ref{_CReadOpts}()
    ccall(_sym(:mio_read_opts_init), Cvoid, (Ptr{_CReadOpts},), base)
    names = opts.arrays === nothing ? String[] : opts.arrays
    roots = [Base.cconvert(Cstring, n) for n in names]
    # One spare slot so that an EMPTY name list still has a non-NULL pointer:
    # NULL means "every array", while a valid pointer with count 0 means "no
    # arrays at all", and that distinction is load-bearing at the ABI.
    ptrs = Cstring[Base.unsafe_convert(Cstring, r) for r in roots]
    push!(ptrs, Cstring(C_NULL))
    GC.@preserve roots ptrs begin
        arrays_ptr = opts.arrays === nothing ? Ptr{Cstring}(C_NULL) : pointer(ptrs)
        cfg = _CReadOpts(opts.points_only ? Cint(1) : Cint(0),
                         opts.metadata_only ? Cint(1) : Cint(0),
                         arrays_ptr, Int64(length(names)), _mmap_mode(opts.mmap),
                         Cint(0), Int64(opts.time_step), base[].reserved)
        ref = Ref(cfg)
        GC.@preserve ref f(ref)
    end
end

"""
    read(path; format="", options=nothing) -> Mesh

Read a mesh file. `format` empty infers from the extension. Pass
[`ReadOptions`](@ref) to narrow what is materialized.

Not exported, because it would shadow `Base.read`: call it qualified, e.g.
`import MeshioPlusPlus as mio; mio.read("bracket.msh")`.
"""
function read(path::AbstractString; format::AbstractString="",
              options::Union{Nothing,ReadOptions}=nothing)
    ptr = if options === nothing
        ccall(_sym(:mio_read), Ptr{Cvoid}, (Cstring, Cstring), path, format)
    else
        _with_read_opts(options) do ref
            ccall(_sym(:mio_read_ex), Ptr{Cvoid}, (Cstring, Cstring, Ptr{_CReadOpts}),
                  path, format, ref)
        end
    end
    Mesh(_check_ptr(ptr))
end

"""
    write(mesh, path; format="")

Write a mesh. `format` empty infers from the extension.

Not exported (it would shadow `Base.write`): call it qualified.
"""
function write(m::Mesh, path::AbstractString; format::AbstractString="")
    _check(ccall(_sym(:mio_write), Cint, (Cstring, Ptr{Cvoid}, Cstring),
                 path, _handle(m), format))
end

"""
    convert(in_path, out_path; in_format="", out_format="")

Read `in_path` and immediately write `out_path`, without materializing a mesh
for the caller (the CLI's `convert`).

Not exported (it would shadow `Base.convert`): call it qualified.
"""
function convert(in_path::AbstractString, out_path::AbstractString;
                 in_format::AbstractString="", out_format::AbstractString="")
    _check(ccall(_sym(:mio_convert), Cint, (Cstring, Cstring, Cstring, Cstring),
                 in_path, in_format, out_path, out_format))
end

"""
    reader_supports_options(format) -> Bool

Whether `format` has a *native* selective-read path. A format without one is
still readable with [`ReadOptions`](@ref) — it is just read in full first.
"""
function reader_supports_options(format::AbstractString)
    r = ccall(_sym(:mio_reader_supports_options), Cint, (Cstring,), format)
    r < 0 && _throw_status(MIO_ERR_INVALID_ARG)
    r != 0
end

"""
    MeshMetadata

Summary of a mesh file produced by [`read_metadata`](@ref) without loading its
heavy arrays. `bbox` is `nothing` for a native summary, which is the normal
case: computing it would mean decoding the coordinates and defeating the
purpose. `fell_back` is `true` when the format had no header-only path, so the
whole file was read (the summary is still correct, just not cheap).
"""
"""One named region's shape, without its entries -- the `read_metadata`
counterpart of [`Region`](@ref). Cheap to enumerate without the cost of
loading every entry."""
struct RegionSummary
    name::String
    kind::Symbol
    dim::Int
    tag::Int
    num_entries::Int
end

struct MeshMetadata
    num_points::Int
    point_dim::Int
    num_cells::Int
    num_cell_blocks::Int
    cell_blocks::Vector{NamedTuple{(:type, :num_cells, :nodes_per_cell, :is_ragged),
                                   Tuple{String,Int,Int,Bool}}}
    point_data_names::Vector{String}
    cell_data_names::Vector{String}
    field_data_names::Vector{String}
    bbox::Union{Nothing,Tuple{NTuple{3,Float64},NTuple{3,Float64}}}
    fell_back::Bool
    """The file's recorded time-series values; empty for a format with no time
    concept. This is the count `ReadOptions(time_step=...)` may name."""
    time_values::Vector{Float64}
    """The file's named regions, without their entries. Empty on a native
    metadata path (none of those formats currently map regions); cheap
    whenever the summary came from an already-read mesh."""
    regions::Vector{RegionSummary}
end

"""
    read_metadata(path; format="") -> MeshMetadata

Summarize a mesh file without loading its heavy arrays.
"""
function read_metadata(path::AbstractString; format::AbstractString="")
    h = _check_ptr(ccall(_sym(:mio_read_metadata_create), Ptr{Cvoid}, (Cstring, Cstring),
                         path, format))
    try
        nblocks = ccall(_sym(:mio_read_metadata_num_cell_blocks), Int64, (Ptr{Cvoid},), h)
        blocks = NamedTuple{(:type, :num_cells, :nodes_per_cell, :is_ragged),
                            Tuple{String,Int,Int,Bool}}[]
        for i in 0:(nblocks-1)
            nc = Ref{Int64}(0); npc = Ref{Int64}(0); ragged = Ref{Cint}(0)
            _check(ccall(_sym(:mio_read_metadata_cell_block), Cint,
                         (Ptr{Cvoid}, Int64, Ptr{Int64}, Ptr{Int64}, Ptr{Cint}),
                         h, i, nc, npc, ragged))
            t = _getstring() do buf, len
                ccall(_sym(:mio_read_metadata_cell_block_type), Int64,
                      (Ptr{Cvoid}, Int64, Ptr{UInt8}, Int64), h, i, buf, len)
            end
            push!(blocks, (type=t, num_cells=Int(nc[]), nodes_per_cell=Int(npc[]),
                           is_ragged=ragged[] != 0))
        end
        names = map((MIO_DATA_POINT, MIO_DATA_CELL, MIO_DATA_FIELD)) do loc
            n = ccall(_sym(:mio_read_metadata_num_names), Int64, (Ptr{Cvoid}, Cint), h, loc)
            n < 0 && (n = 0)
            String[_getstring() do buf, len
                       ccall(_sym(:mio_read_metadata_name), Int64,
                             (Ptr{Cvoid}, Cint, Int64, Ptr{UInt8}, Int64), h, loc, i, buf, len)
                   end for i in 0:(n-1)]
        end
        lo = Vector{Float64}(undef, 3); hi = Vector{Float64}(undef, 3)
        st = ccall(_sym(:mio_read_metadata_bbox), Cint,
                   (Ptr{Cvoid}, Ptr{Float64}, Ptr{Float64}), h, lo, hi)
        bbox = st == MIO_OK ? ((lo[1], lo[2], lo[3]), (hi[1], hi[2], hi[3])) : nothing
        nsteps = ccall(_sym(:mio_read_metadata_num_time_values), Int64, (Ptr{Cvoid},), h)
        nsteps < 0 && (nsteps = 0)
        times = Vector{Float64}(undef, nsteps)
        nsteps > 0 && ccall(_sym(:mio_read_metadata_time_values), Int64,
                            (Ptr{Cvoid}, Ptr{Float64}, Int64), h, times, nsteps)
        nregions = ccall(_sym(:mio_read_metadata_num_regions), Int64, (Ptr{Cvoid},), h)
        nregions < 0 && (nregions = 0)
        rsummaries = RegionSummary[]
        for i in 0:(nregions-1)
            rname = _getstring() do buf, len
                ccall(_sym(:mio_read_metadata_region_name), Int64,
                      (Ptr{Cvoid}, Int64, Ptr{UInt8}, Int64), h, i, buf, len)
            end
            info = Ref{_CRegionInfo}()
            _check(ccall(_sym(:mio_read_metadata_region_info), Cint,
                         (Ptr{Cvoid}, Int64, Ptr{_CRegionInfo}), h, i, info))
            ci = info[]
            push!(rsummaries,
                  RegionSummary(rname, _kind_sym(ci.kind), Int(ci.dim), Int(ci.tag),
                               Int(ci.num_entries)))
        end
        MeshMetadata(
            Int(ccall(_sym(:mio_read_metadata_num_points), Int64, (Ptr{Cvoid},), h)),
            Int(ccall(_sym(:mio_read_metadata_point_dim), Int64, (Ptr{Cvoid},), h)),
            Int(ccall(_sym(:mio_read_metadata_num_cells), Int64, (Ptr{Cvoid},), h)),
            Int(nblocks), blocks, names[1], names[2], names[3], bbox,
            ccall(_sym(:mio_read_metadata_fell_back), Cint, (Ptr{Cvoid},), h) == 1,
            times, rsummaries)
    finally
        ccall(_sym(:mio_read_metadata_free), Cvoid, (Ptr{Cvoid},), h)
    end
end

# --- inspecting: counts ------------------------------------------------------

"""    num_points(mesh) -> Int"""
num_points(m::Mesh) = Int(ccall(_sym(:mio_mesh_num_points), Int64, (Ptr{Cvoid},), _handle(m)))

"""    point_dim(mesh) -> Int — coordinates per point (usually 2 or 3)."""
point_dim(m::Mesh) = Int(ccall(_sym(:mio_mesh_point_dim), Int64, (Ptr{Cvoid},), _handle(m)))

"""    num_cell_blocks(mesh) -> Int"""
num_cell_blocks(m::Mesh) =
    Int(ccall(_sym(:mio_mesh_num_cell_blocks), Int64, (Ptr{Cvoid},), _handle(m)))

"""
    num_cells(mesh) -> Int

Total cells over all blocks. The C API has no single accessor for this (a
documented gap), so it is summed over [`cell_block_info`](@ref) here.
"""
function num_cells(m::Mesh)
    total = 0
    for i in 1:num_cell_blocks(m)
        total += cell_block_info(m, i).num_cells
    end
    total
end

"""
    cell_block_info(mesh, block) -> (; num_cells, nodes_per_cell, is_ragged)

Shape of cell block `block` (**1-based**). A ragged block (polygons or
polyhedra of varying size) reports `nodes_per_cell == 0` and `is_ragged ==
true`; its connectivity is **not reachable through the C API** (a documented
v1 limitation), so [`connectivity`](@ref) throws for it.
"""
function cell_block_info(m::Mesh, block::Integer)
    nc = Ref{Int64}(0); npc = Ref{Int64}(0); ragged = Ref{Int32}(0)
    _check(ccall(_sym(:mio_mesh_cell_block_info), Cint,
                 (Ptr{Cvoid}, Int64, Ptr{Int64}, Ptr{Int64}, Ptr{Int32}),
                 _handle(m), Int64(block - 1), nc, npc, ragged))
    (num_cells=Int(nc[]), nodes_per_cell=Int(npc[]), is_ragged=ragged[] != 0)
end

"""    cell_block_type(mesh, block) -> String — the meshio type name (1-based `block`)."""
function cell_block_type(m::Mesh, block::Integer)
    h = _handle(m)
    _getstring() do buf, len
        ccall(_sym(:mio_mesh_cell_block_type), Int64,
              (Ptr{Cvoid}, Int64, Ptr{UInt8}, Int64), h, Int64(block - 1), buf, len)
    end
end

"""    cell_block_types(mesh) -> Vector{String}"""
cell_block_types(m::Mesh) = [cell_block_type(m, i) for i in 1:num_cell_blocks(m)]

# --- inspecting: points ------------------------------------------------------

"""
    points_ptr(mesh) -> MeshBorrow

**Zero-copy** `(dim, num_points)` view of the coordinates, in whatever dtype
the mesh stores (`Float32` or `Float64`). Column-major, so this is the same
memory as the C API's row-major `(num_points, dim)` — nothing is transposed.

Valid until the next mutating call on `mesh` or [`close`](@ref); see
[`MeshBorrow`](@ref). Coordinates carry no indices, so there is no 0-/1-based
question here.
"""
function points_ptr(m::Mesh)
    h = _handle(m)
    data = Ref{Ptr{Cvoid}}(C_NULL); dt = Ref{Cint}(0)
    _check(ccall(_sym(:mio_mesh_get_points), Cint,
                 (Ptr{Cvoid}, Ptr{Ptr{Cvoid}}, Ptr{Cint}), h, data, dt))
    _borrow(m, data[], _jltype(dt[]), (point_dim(m), num_points(m)))
end

"""
    points(mesh) -> Matrix

**Copy** of the `(dim, num_points)` coordinates, in the mesh's own dtype. Use
[`points_ptr`](@ref) to avoid the copy.
"""
points(m::Mesh) = copy(parent(points_ptr(m)))

# --- inspecting: connectivity ------------------------------------------------

"""
    connectivity_ptr(mesh, block) -> MeshBorrow

**Zero-copy, 0-BASED** `(nodes_per_cell, num_cells)` view of cell block
`block`'s connectivity (`block` itself is 1-based, like every index in this
package). The node indices are the ABI's own 0-based ones: a borrow cannot be
shifted without copying it, which is exactly what a borrow exists to avoid.

Use [`connectivity`](@ref) for a 1-based copy. Throws on a ragged block, whose
connectivity the C API does not expose (v1).
"""
function connectivity_ptr(m::Mesh, block::Integer)
    h = _handle(m)
    info = cell_block_info(m, block)
    conn = Ref{Ptr{Cvoid}}(C_NULL); dt = Ref{Cint}(0)
    _check(ccall(_sym(:mio_mesh_cell_block_conn), Cint,
                 (Ptr{Cvoid}, Int64, Ptr{Ptr{Cvoid}}, Ptr{Cint}),
                 h, Int64(block - 1), conn, dt))
    _borrow(m, conn[], _jltype(dt[]), (info.nodes_per_cell, info.num_cells))
end

"""
    connectivity(mesh, block) -> Matrix{Int64}

**Copy, 1-BASED** `(nodes_per_cell, num_cells)` connectivity of cell block
`block` (1-based). The `+1` shift happens here, where a copy is made anyway;
[`connectivity_ptr`](@ref) is the un-shifted zero-copy borrow.

Throws on a ragged block: the C API reports `is_ragged` but does not expose
that connectivity (a documented v1 limitation).
"""
function connectivity(m::Mesh, block::Integer)
    b = connectivity_ptr(m, block)
    out = Matrix{Int64}(undef, size(b))
    src = parent(b)
    @inbounds for i in eachindex(src, out)
        out[i] = Int64(src[i]) + 1
    end
    out
end

# --- inspecting: data arrays -------------------------------------------------

# The C API reports a data array as row-major `ndim`/`shape`; the column-major
# Julia array over the same memory has the reversed shape. A scalar per-point
# field is a length-n vector either way; an (n, 3) vector field is (3, n) here.
function _data_dims(ndim::Int32, shape::Vector{Int64})
    Tuple(reverse(shape[1:ndim]))
end

for (kind, cnum, cname) in ((:point_data, :mio_mesh_num_point_data, :mio_mesh_point_data_name),
                            (:cell_data, :mio_mesh_num_cell_data, :mio_mesh_cell_data_name),
                            (:field_data, :mio_mesh_num_field_data, :mio_mesh_field_data_name))
    names_fn = Symbol(kind, "_names")
    count_fn = Symbol("num_", kind)
    @eval begin
        """    $($(string(count_fn)))(mesh) -> Int"""
        $count_fn(m::Mesh) = Int(ccall(_sym($(QuoteNode(cnum))), Int64, (Ptr{Cvoid},), _handle(m)))

        """
            $($(string(names_fn)))(mesh) -> Vector{String}

        The array names, in ascending lexicographic order — identical on every
        mesh backend, which is what keeps output byte-identical across them.
        """
        function $names_fn(m::Mesh)
            h = _handle(m)
            n = $count_fn(m)
            String[_getstring() do buf, len
                       ccall(_sym($(QuoteNode(cname))), Int64,
                             (Ptr{Cvoid}, Int64, Ptr{UInt8}, Int64), h, Int64(i), buf, len)
                   end for i in 0:(n-1)]
        end
    end
end

"""
    point_data_ptr(mesh, name) -> MeshBorrow

**Zero-copy** view of the named per-point array, shaped `(components...,
num_points)` — the column-major mirror of the C API's row-major
`(num_points, components...)`.
"""
function point_data_ptr(m::Mesh, name::AbstractString)
    h = _handle(m)
    data = Ref{Ptr{Cvoid}}(C_NULL); dt = Ref{Cint}(0); nd = Ref{Int32}(0)
    shape = zeros(Int64, MIO_MAX_NDIM)
    _check(ccall(_sym(:mio_mesh_get_point_data), Cint,
                 (Ptr{Cvoid}, Cstring, Ptr{Ptr{Cvoid}}, Ptr{Cint}, Ptr{Int32}, Ptr{Int64}),
                 h, name, data, dt, nd, shape))
    _borrow(m, data[], _jltype(dt[]), _data_dims(nd[], shape))
end

"""    point_data(mesh, name) -> Array — copy of [`point_data_ptr`](@ref)."""
point_data(m::Mesh, name::AbstractString) = copy(parent(point_data_ptr(m, name)))

"""
    cell_data_num_blocks(mesh, name) -> Int

How many per-block arrays the named cell-data field has (normally the number
of cell blocks).
"""
function cell_data_num_blocks(m::Mesh, name::AbstractString)
    n = ccall(_sym(:mio_mesh_cell_data_num_blocks), Int64, (Ptr{Cvoid}, Cstring),
              _handle(m), name)
    Int(_check_count(n, "cell_data_num_blocks(\"$name\")"))
end

"""
    cell_data_ptr(mesh, name, block) -> MeshBorrow

**Zero-copy** view of the named cell-data field's array for cell block
`block` (1-based), shaped `(components..., num_cells_in_block)`.
"""
function cell_data_ptr(m::Mesh, name::AbstractString, block::Integer)
    h = _handle(m)
    data = Ref{Ptr{Cvoid}}(C_NULL); dt = Ref{Cint}(0); nd = Ref{Int32}(0)
    shape = zeros(Int64, MIO_MAX_NDIM)
    _check(ccall(_sym(:mio_mesh_get_cell_data), Cint,
                 (Ptr{Cvoid}, Cstring, Int64, Ptr{Ptr{Cvoid}}, Ptr{Cint}, Ptr{Int32},
                  Ptr{Int64}), h, name, Int64(block - 1), data, dt, nd, shape))
    _borrow(m, data[], _jltype(dt[]), _data_dims(nd[], shape))
end

"""    cell_data(mesh, name, block) -> Array — copy of [`cell_data_ptr`](@ref)."""
cell_data(m::Mesh, name::AbstractString, block::Integer) =
    copy(parent(cell_data_ptr(m, name, block)))

"""
    field_data_ptr(mesh, name) -> MeshBorrow

**Zero-copy** view of the named mesh-level (field) array, with the C shape
reversed for column-major order.
"""
function field_data_ptr(m::Mesh, name::AbstractString)
    h = _handle(m)
    data = Ref{Ptr{Cvoid}}(C_NULL); dt = Ref{Cint}(0); nd = Ref{Int32}(0)
    shape = zeros(Int64, MIO_MAX_NDIM)
    _check(ccall(_sym(:mio_mesh_get_field_data), Cint,
                 (Ptr{Cvoid}, Cstring, Ptr{Ptr{Cvoid}}, Ptr{Cint}, Ptr{Int32}, Ptr{Int64}),
                 h, name, data, dt, nd, shape))
    _borrow(m, data[], _jltype(dt[]), _data_dims(nd[], shape))
end

"""    field_data(mesh, name) -> Array — copy of [`field_data_ptr`](@ref)."""
field_data(m::Mesh, name::AbstractString) = copy(parent(field_data_ptr(m, name)))

# --- building (setters -- all COPY caller memory) ----------------------------

# Dense column-major array of element type `T`, copying only when needed. A
# plain `Array` is always dense, so it is passed straight through; anything
# else (a view, an adjoint, a different element type) is materialized first,
# because the C API reads the buffer as a flat row-major block.
_dense(A::Array{T}, ::Type{T}) where {T} = A
_dense(A::AbstractArray, ::Type{T}) where {T} = Array{T}(A)

"""
    set_points!(mesh, points)

Assign the coordinates from a `(dim, num_points)` matrix, replacing any
previous ones. `Float32` and `Float64` are stored as given.

Column-major: the matrix is the same memory as the C API's row-major
`(num_points, dim)`, so nothing is transposed. **Copies** (C API rule 2) and
invalidates every outstanding borrow into `mesh`.
"""
function set_points!(m::Mesh, pts::AbstractMatrix{T}) where {T<:AbstractFloat}
    T in (Float32, Float64) ||
        throw(ArgumentError("points must be Float32 or Float64, got $T"))
    A = _dense(pts, T)
    dim, n = size(A)
    _check(GC.@preserve A ccall(_sym(:mio_mesh_set_points), Cint,
                                (Ptr{Cvoid}, Cint, Int64, Int64, Ptr{Cvoid}),
                                _handle(m), _miodtype(T), Int64(n), Int64(dim), pointer(A)))
    _touch!(m)
    m
end

"""
    add_cell_block!(mesh, cell_type, conn)

Append one homogeneous cell block from a **1-based** `(nodes_per_cell,
num_cells)` connectivity matrix. The `-1` shift to the core's 0-based indexing
happens here, in the copy.

`cell_type` is a meshio type name (`"triangle"`, `"tetra10"`, …). **Copies**
and invalidates every outstanding borrow into `mesh`.
"""
function add_cell_block!(m::Mesh, cell_type::AbstractString,
                         conn::AbstractMatrix{<:Integer})
    npc, nc = size(conn)
    shifted = Matrix{Int64}(undef, npc, nc)
    @inbounds for i in eachindex(conn, shifted)
        v = Int64(conn[i])
        v >= 1 || throw(ArgumentError(
            "connectivity is 1-based here; got $v (use the 0-based ABI only " *
            "through connectivity_ptr)"))
        shifted[i] = v - 1
    end
    _check(GC.@preserve shifted ccall(_sym(:mio_mesh_add_cell_block), Cint,
                                      (Ptr{Cvoid}, Cstring, Int64, Int64, Cint, Ptr{Cvoid}),
                                      _handle(m), cell_type, Int64(nc), Int64(npc),
                                      _miodtype(Int64), pointer(shifted)))
    _touch!(m)
    m
end

# Shared marshalling for the three data setters: a Julia array shaped
# (components..., rows) is the C row-major (rows, components...).
function _add_data(m::Mesh, sym::Symbol, name::AbstractString,
                   data::AbstractArray{T}) where {T}
    A = _dense(data, T)
    shape = Int64[reverse(size(A))...]
    ndim = Int32(length(shape))
    st = GC.@preserve A shape ccall(_sym(sym), Cint,
                                    (Ptr{Cvoid}, Cstring, Cint, Int32, Ptr{Int64}, Ptr{Cvoid}),
                                    _handle(m), name, _miodtype(T), ndim, pointer(shape),
                                    pointer(A))
    _check(st)
    _touch!(m)
    m
end

"""
    add_point_data!(mesh, name, data)

Attach a named per-point array. `data` is `(components..., num_points)`: a
length-`num_points` vector for a scalar field, a `(3, num_points)` matrix for a
vector field. **Copies** and invalidates outstanding borrows.
"""
add_point_data!(m::Mesh, name::AbstractString, data::AbstractArray) =
    _add_data(m, :mio_mesh_add_point_data, name, data)

"""
    append_cell_data!(mesh, name, data)

Append one per-cell-block array to the named cell-data field. Call it **once
per cell block, in block order**, after adding the blocks; the trailing
dimension must equal that block's cell count. **Copies** and invalidates
outstanding borrows.
"""
append_cell_data!(m::Mesh, name::AbstractString, data::AbstractArray) =
    _add_data(m, :mio_mesh_append_cell_data, name, data)

"""
    add_field_data!(mesh, name, data)

Attach a named mesh-level (field) array of arbitrary shape. **Copies** and
invalidates outstanding borrows.
"""
add_field_data!(m::Mesh, name::AbstractString, data::AbstractArray) =
    _add_data(m, :mio_mesh_add_field_data, name, data)
