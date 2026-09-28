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

const _REGION_OPS = Dict(:union => Int32(0), :intersection => Int32(1),
                         :intersect => Int32(1), :difference => Int32(2),
                         :rename => Int32(3), :retag => Int32(4), :delete => Int32(5))

"""
    edit_regions(mesh, op, inputs; output="", kind=nothing, dim=nothing, tag=nothing,
                 keep_inputs=true) -> Mesh

Apply one region edit to a copy of `mesh` (points, cells and data untouched).
`op` is `:union`, `:intersection`, `:difference` (two or more `inputs` of one
kind, the result named `output`), `:rename` (one input, to `output`), `:retag`
(one input, new `tag` and/or `dim`) or `:delete`. `inputs` are region names;
`kind` (`:point`, `:cell`, `:side`) pins their kind when a name is shared. See
`doc/regions.md`.
"""
function edit_regions(m::Mesh, op::Symbol, inputs::AbstractVector{<:AbstractString};
                      output::AbstractString="", kind::Union{Nothing,Symbol}=nothing,
                      dim::Union{Nothing,Integer}=nothing,
                      tag::Union{Nothing,Integer}=nothing, keep_inputs::Bool=true)
    haskey(_REGION_OPS, op) ||
        throw(ArgumentError("meshio++: edit_regions: unknown operation :$op"))
    names = [Vector{UInt8}(codeunits(String(n) * "\0")) for n in inputs]
    out_c = Vector{UInt8}(codeunits(String(output) * "\0"))
    k = kind === nothing ? Int32(-1) : Int32(_kind_val(kind))
    ptr = GC.@preserve names out_c begin
        sels = [_CRegionSelector(Cstring(pointer(n)), k, Int32(0), Int64(-2), Int64(-2),
                                 (Int64(0), Int64(0))) for n in names]
        ccall(_sym(:mio_edit_regions), Ptr{Cvoid},
              (Ptr{Cvoid}, Int32, Ptr{_CRegionSelector}, Int64, Cstring, Int64, Int64, Int32),
              _handle(m), _REGION_OPS[op], sels, Int64(length(sels)),
              isempty(output) ? C_NULL : Cstring(pointer(out_c)),
              dim === nothing ? Int64(-2) : Int64(dim),
              tag === nothing ? Int64(-2) : Int64(tag), Int32(keep_inputs))
    end
    Mesh(_check_ptr(ptr))
end

"""
    remove_region!(mesh, index)

Remove the `index`-th region (1-based, [`regions`](@ref) order). Invalidates
every outstanding borrow into `mesh`.
"""
function remove_region!(m::Mesh, index::Integer)
    _check(ccall(_sym(:mio_mesh_remove_region), Cint, (Ptr{Cvoid}, Int64),
                 _handle(m), Int64(index - 1)))
    m
end

"""
    match_periodic_nodes(mesh, slave, master; translate=nothing, matrix=nothing,
                         atol=1e-8, require_complete=true)
        -> (; slave, master, unmatched, num_fixed, max_residual)

The master node each node of the `slave` region maps onto under an affine
transform -- `translate` (3 numbers) or `matrix` (16 numbers, a row-major 4x4)
-- within `atol`. Node ids are **1-based**; `slave` is ascending and `master`
aligned with it. With `require_complete=false` unmatched slave nodes are
returned instead of failing. See `doc/periodic.md`.
"""
function match_periodic_nodes(m::Mesh, slave::AbstractString, master::AbstractString;
                              translate=nothing, matrix=nothing, atol::Real=1e-8,
                              require_complete::Bool=true)
    mat = zeros(16)
    mat[1] = mat[6] = mat[11] = mat[16] = 1.0
    if matrix !== nothing
        length(matrix) == 16 ||
            throw(ArgumentError("meshio++: match_periodic_nodes: matrix needs 16 numbers"))
        mat .= Float64.(vec(collect(matrix)))
    elseif translate !== nothing
        length(translate) == 3 ||
            throw(ArgumentError("meshio++: match_periodic_nodes: translate needs 3 numbers"))
        mat[4], mat[8], mat[12] = Float64.(collect(translate))
    end
    s_c = Vector{UInt8}(codeunits(String(slave) * "\0"))
    m_c = Vector{UInt8}(codeunits(String(master) * "\0"))
    handle = GC.@preserve s_c m_c begin
        ss = _CRegionSelector(Cstring(pointer(s_c)), Int32(-1), Int32(0), Int64(-2), Int64(-2),
                              (Int64(0), Int64(0)))
        ms = _CRegionSelector(Cstring(pointer(m_c)), Int32(-1), Int32(0), Int64(-2), Int64(-2),
                              (Int64(0), Int64(0)))
        opts = _CPeriodicOpts(Tuple(mat), Float64(atol), Int32(require_complete), Int32(0),
                              ntuple(_ -> Int64(0), 6))
        ccall(_sym(:mio_match_periodic_nodes), Ptr{Cvoid},
              (Ptr{Cvoid}, Ref{_CRegionSelector}, Ref{_CRegionSelector}, Ref{_CPeriodicOpts}),
              _handle(m), Ref(ss), Ref(ms), Ref(opts))
    end
    _check_ptr(handle)
    try
        np, nu, nf = Ref{Int64}(0), Ref{Int64}(0), Ref{Int64}(0)
        res = Ref{Cdouble}(0.0)
        _check(ccall(_sym(:mio_periodic_pairs_info), Cint,
                     (Ptr{Cvoid}, Ptr{Int64}, Ptr{Int64}, Ptr{Int64}, Ptr{Cdouble}),
                     handle, np, nu, nf, res))
        take(sym, n) = begin
            n == 0 && return Int64[]
            cnt = Ref{Int64}(0)
            p = ccall(_sym(sym), Ptr{Int64}, (Ptr{Cvoid}, Ptr{Int64}), handle, cnt)
            p == C_NULL && _throw_status(MIO_ERR_INTERNAL)
            unsafe_wrap(Array, p, (Int(cnt[]),); own=false) .+ 1
        end
        (slave=take(:mio_periodic_pairs_slave, np[]),
         master=take(:mio_periodic_pairs_master, np[]),
         unmatched=take(:mio_periodic_pairs_unmatched, nu[]),
         num_fixed=Int(nf[]), max_residual=res[])
    finally
        ccall(_sym(:mio_periodic_pairs_free), Cvoid, (Ptr{Cvoid},), handle)
    end
end
