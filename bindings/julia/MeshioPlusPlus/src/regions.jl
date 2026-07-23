# Named regions -- the unified model of a named group of mesh entities (a gmsh
# physical group, an Exodus block / node set / side set, an Abaqus *NSET /
# *ELSET / *SURFACE, a MED family, a Kratos SubModelPart). See doc/regions.md.

"""
    Region(name, kind, dim, tag, entries)

A named group of mesh entities.

* `kind` is `:point` (point indices), `:cell` (**global, block-major** cell
  indices — block 1's cells first, then block 2's, …) or `:side`
  (`(cell, facet)` pairs).
* `dim` is the topological dimension, or `-1` when unspecified; `tag` is the
  format-native integer id, or `-1` when there is none. Both exist because a
  gmsh physical group is per-dimension and carries a tag that two groups of
  different dimensions may share.
* `entries` is a `(stride, num_entries)` matrix — `stride` is 2 for `:side`
  and 1 otherwise. Column-major, so it is the same layout as the C API's
  row-major entry buffer.

**Indexing:** point and cell indices are **1-based** here, shifted in the copy
like all the other copying accessors. For `:side`, only row 1 (the cell index)
is shifted — **row 2, the facet number, is left alone**, because a facet is an
ordinal within a cell type, not an index into the mesh. That is the same rule
the Fortran module follows.

Entries come back sorted and de-duplicated: the core canonicalizes every
region on insertion, which is what makes region equality exact and output
byte-identical across mesh backends and thread counts.
"""
struct Region
    name::String
    kind::Symbol
    dim::Int
    tag::Int
    entries::Matrix{Int64}
end

function Base.show(io::IO, r::Region)
    print(io, "Region(\"", r.name, "\", :", r.kind, ", dim=", r.dim, ", tag=", r.tag,
          ", ", size(r.entries, 2), " entries)")
end

const _KIND_SYMS = (:point, :cell, :side)

function _kind_sym(k::Integer)
    0 <= k < length(_KIND_SYMS) ||
        throw(MeshioError(MIO_ERR_INVALID_ARG, "unknown region kind $k"))
    _KIND_SYMS[k+1]
end

function _kind_val(k::Symbol)
    i = findfirst(==(k), _KIND_SYMS)
    i === nothing && throw(ArgumentError(
        "region kind must be :point, :cell or :side, got :$k"))
    Cint(i - 1)
end

"""
    regions(mesh) -> Vector{Region}

Every named region the mesh carries, as **copies** with 1-based indices (see
[`Region`](@ref) for the `:side` facet exception). The C API hands the entries
back as a borrow that would die at the next mutating call, so they are copied
here rather than exposed as a [`MeshBorrow`](@ref).
"""
function regions(m::Mesh)
    h = _handle(m)
    handle = _check_ptr(ccall(_sym(:mio_regions_create), Ptr{Cvoid}, (Ptr{Cvoid},), h))
    try
        n = _check_count(ccall(_sym(:mio_regions_count), Int64, (Ptr{Cvoid},), handle),
                         "mio_regions_count")
        out = Region[]
        for i in 0:(n-1)
            name = _getstring() do buf, len
                ccall(_sym(:mio_regions_name), Int64,
                      (Ptr{Cvoid}, Int64, Ptr{UInt8}, Int64), handle, i, buf, len)
            end
            info = Ref{_CRegionInfo}()
            _check(ccall(_sym(:mio_regions_info), Cint,
                         (Ptr{Cvoid}, Int64, Ptr{_CRegionInfo}), handle, i, info))
            ci = info[]
            kind = _kind_sym(ci.kind)
            stride = Int(ci.stride)
            nent = Int(ci.num_entries)
            entries = Matrix{Int64}(undef, stride, nent)
            if nent > 0
                count = Ref{Int64}(0)
                ptr = ccall(_sym(:mio_regions_entries), Ptr{Int64},
                            (Ptr{Cvoid}, Int64, Ptr{Int64}), handle, i, count)
                ptr == C_NULL && _throw_status(MIO_ERR_INTERNAL)
                src = unsafe_wrap(Array, ptr, (stride, nent); own=false)
                # +1 on the index rows; for :side the facet row (2) is an
                # ordinal within the cell type, not a mesh index, so it stays.
                @inbounds for c in 1:nent, r in 1:stride
                    entries[r, c] = (kind === :side && r == 2) ? src[r, c] : src[r, c] + 1
                end
            end
            push!(out, Region(name, kind, Int(ci.dim), Int(ci.tag), entries))
        end
        out
    finally
        ccall(_sym(:mio_regions_free), Cvoid, (Ptr{Cvoid},), handle)
    end
end

"""
    add_region!(mesh, name, kind, entries; dim=-1, tag=-1)

Add a named region, replacing any existing one with the same
`(kind, name, dim, tag)`.

`entries` is either a `(stride, num_entries)` matrix or, for `:point`/`:cell`,
a plain vector of indices. Indices are **1-based** and shifted back in the
copy; for `:side` the facet row is passed through unshifted (see
[`Region`](@ref)). The core canonicalizes the entries — sorts and
de-duplicates them — so they come back out ordered regardless of the order
given here.

**Copies** and invalidates every outstanding borrow into `mesh`.
"""
function add_region!(m::Mesh, name::AbstractString, kind::Symbol,
                     entries::AbstractMatrix{<:Integer}; dim::Integer=-1,
                     tag::Integer=-1)
    kv = _kind_val(kind)
    stride, nent = size(entries)
    expected = kind === :side ? 2 : 1
    stride == expected || throw(ArgumentError(
        "a :$kind region needs $expected value(s) per entry, got $stride"))
    flat = Vector{Int64}(undef, stride * nent)
    k = 0
    @inbounds for c in 1:nent, r in 1:stride
        v = Int64(entries[r, c])
        if kind === :side && r == 2
            flat[k+=1] = v            # facet ordinal: not a mesh index
        else
            v >= 1 || throw(ArgumentError(
                "region entries are 1-based here; got $v"))
            flat[k+=1] = v - 1
        end
    end
    _check(GC.@preserve flat ccall(_sym(:mio_mesh_add_region), Cint,
                                   (Ptr{Cvoid}, Cstring, Cint, Int32, Int64, Ptr{Int64},
                                    Int64), _handle(m), name, kv, Int32(dim), Int64(tag),
                                   pointer(flat), Int64(length(flat))))
    _touch!(m)
    m
end

function add_region!(m::Mesh, name::AbstractString, kind::Symbol,
                     entries::AbstractVector{<:Integer}; kwargs...)
    kind === :side && throw(ArgumentError(
        "a :side region needs (cell, facet) pairs: pass a 2xN matrix"))
    add_region!(m, name, kind, reshape(collect(entries), 1, :); kwargs...)
end
