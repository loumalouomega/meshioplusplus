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
    region_adjacency(mesh; regions=nothing) -> Mesh

Return the conforming facets shared by selected Cell regions. With no names,
all Cell regions are used; if fewer than two exist, cell blocks become groups.
The result carries parent cell/facet ids and per-facet length or area. See
`doc/region_adjacency.md`.
"""
function region_adjacency(m::Mesh; regions::Union{Nothing,AbstractVector{<:AbstractString}}=nothing)
    names = regions === nothing ? String[] : String.(collect(regions))
    names_c = [Vector{UInt8}(codeunits(name * "\0")) for name in names]
    ptr = GC.@preserve names_c begin
        selectors = [_CRegionSelector(Cstring(pointer(name)), Int32(1), Int32(0),
                                      Int64(-2), Int64(-2), (Int64(0), Int64(0)))
                     for name in names_c]
        ccall(_sym(:mio_region_adjacency), Ptr{Cvoid},
              (Ptr{Cvoid}, Ptr{_CRegionSelector}, Int64),
              _handle(m), selectors, Int64(length(selectors)))
    end
    Mesh(_check_ptr(ptr))
end

"""
    find_interface(mesh_a, region_a, region_b; mesh_b=nothing, mode=:conforming,
                   master=:a, gap_tolerance=0, angle_tolerance=30,
                   overlap_tolerance=0)

Find conforming or proximity-matched interfaces. Returns a named tuple with the
master-side facet `mesh` and a `report` containing counts, measure and 0-based
`(cell, local_facet)` Side entries for each source mesh. Cell ids in Julia are
1-based; local facet ordinals remain 0-based. See `doc/region_adjacency.md`.
"""
function find_interface(a::Mesh, region_a::AbstractString, region_b::AbstractString;
                        mesh_b::Union{Nothing,Mesh}=nothing, mode::Symbol=:conforming,
                        master::Symbol=:a, gap_tolerance::Real=0,
                        angle_tolerance::Real=30, overlap_tolerance::Real=0)
    mode_code = mode === :conforming ? Int32(0) : mode === :proximity ? Int32(1) :
        throw(ArgumentError("mode must be :conforming or :proximity"))
    master_code = master === :a ? Int32(0) : master === :b ? Int32(1) :
        throw(ArgumentError("master must be :a or :b"))
    a_c = Vector{UInt8}(codeunits(String(region_a) * "\0"))
    b_c = Vector{UInt8}(codeunits(String(region_b) * "\0"))
    handle = GC.@preserve a_c b_c begin
        sa = _CRegionSelector(Cstring(pointer(a_c)), Int32(1), Int32(0), Int64(-2),
                              Int64(-2), (Int64(0), Int64(0)))
        sb = _CRegionSelector(Cstring(pointer(b_c)), Int32(1), Int32(0), Int64(-2),
                              Int64(-2), (Int64(0), Int64(0)))
        opts = _CFindInterfaceOpts(mode_code, master_code, Float64(gap_tolerance),
                                   Float64(angle_tolerance), Float64(overlap_tolerance),
                                   ntuple(_ -> Int64(0), 4))
        ccall(_sym(:mio_find_interface), Ptr{Cvoid},
              (Ptr{Cvoid}, Ref{_CRegionSelector}, Ptr{Cvoid}, Ref{_CRegionSelector},
               Ref{_CFindInterfaceOpts}), _handle(a), Ref(sa),
              mesh_b === nothing ? C_NULL : _handle(mesh_b), Ref(sb), Ref(opts))
    end
    _check_ptr(handle)
    try
        pairs, unmatched_a, unmatched_b = Ref{Int64}(0), Ref{Int64}(0), Ref{Int64}(0)
        area, max_gap = Ref{Cdouble}(0), Ref{Cdouble}(0)
        _check(ccall(_sym(:mio_find_interface_result_report), Cint,
                     (Ptr{Cvoid}, Ptr{Int64}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Int64}, Ptr{Int64}),
                     handle, pairs, area, max_gap, unmatched_a, unmatched_b))
        function side(sym)
            count = Ref{Int64}(0)
            ptr = ccall(_sym(sym), Ptr{Int64}, (Ptr{Cvoid}, Ptr{Int64}), handle, count)
            n = Int(count[])
            n == 0 && return zeros(Int64, 2, 0)
            ptr == C_NULL && _throw_status(MIO_ERR_INTERNAL)
            entries = copy(unsafe_wrap(Array, ptr, (2, n); own=false))
            entries[1, :] .+= 1 # global cell ids are 1-based in Julia
            entries
        end
        mesh_ptr = ccall(_sym(:mio_find_interface_result_take_mesh), Ptr{Cvoid},
                         (Ptr{Cvoid},), handle)
        mesh = Mesh(_check_ptr(mesh_ptr))
        report = (num_pairs=Int(pairs[]), area=area[], max_gap=max_gap[],
                  unmatched_a=Int(unmatched_a[]), unmatched_b=Int(unmatched_b[]),
                  side_a=side(:mio_find_interface_result_side_a),
                  side_b=side(:mio_find_interface_result_side_b))
        (mesh=mesh, report=report)
    finally
        ccall(_sym(:mio_find_interface_result_free), Cvoid, (Ptr{Cvoid},), handle)
    end
end

"""
    contact_pairs(slave_mesh, slave_points, master_cells; master_mesh=nothing,
                  tolerance=0, require_complete=false)

Project a Point region to the closest facets of a Cell region. Point and global
cell ids are returned 1-based; local facet/subfacet ordinals remain 0-based.
Array results are copied into Julia-owned storage. See `doc/region_adjacency.md`.
"""
function contact_pairs(slave::Mesh, slave_points::AbstractString, master_cells::AbstractString;
                       master_mesh::Union{Nothing,Mesh}=nothing, tolerance::Real=0,
                       require_complete::Bool=false)
    slave_c = Vector{UInt8}(codeunits(String(slave_points) * "\0"))
    master_c = Vector{UInt8}(codeunits(String(master_cells) * "\0"))
    handle = GC.@preserve slave_c master_c begin
        ss = _CRegionSelector(Cstring(pointer(slave_c)), Int32(0), Int32(0), Int64(-2),
                              Int64(-2), (Int64(0), Int64(0)))
        ms = _CRegionSelector(Cstring(pointer(master_c)), Int32(1), Int32(0), Int64(-2),
                              Int64(-2), (Int64(0), Int64(0)))
        opts = _CContactPairsOpts(Float64(tolerance), Int32(require_complete), Int32(0),
                                  ntuple(_ -> Int64(0), 4))
        ccall(_sym(:mio_contact_pairs), Ptr{Cvoid},
              (Ptr{Cvoid}, Ref{_CRegionSelector}, Ptr{Cvoid}, Ref{_CRegionSelector},
               Ref{_CContactPairsOpts}), _handle(slave), Ref(ss),
              master_mesh === nothing ? C_NULL : _handle(master_mesh), Ref(ms), Ref(opts))
    end
    _check_ptr(handle)
    try
        count, unmatched_count = Ref{Int64}(0), Ref{Int64}(0)
        _check(ccall(_sym(:mio_contact_pairs_result_info), Cint,
                     (Ptr{Cvoid}, Ptr{Int64}, Ptr{Int64}), handle, count, unmatched_count))
        n, nu = Int(count[]), Int(unmatched_count[])
        take_i64(sym, len) = len == 0 ? Int64[] : begin
            ptr = ccall(_sym(sym), Ptr{Int64}, (Ptr{Cvoid},), handle)
            ptr == C_NULL && _throw_status(MIO_ERR_INTERNAL)
            copy(unsafe_wrap(Array, ptr, (len,); own=false))
        end
        take_f64(sym, len) = len == 0 ? Float64[] : begin
            ptr = ccall(_sym(sym), Ptr{Cdouble}, (Ptr{Cvoid},), handle)
            ptr == C_NULL && _throw_status(MIO_ERR_INTERNAL)
            copy(unsafe_wrap(Array, ptr, (len,); own=false))
        end
        slave_ids = take_i64(:mio_contact_pairs_slave_point, n) .+ 1
        master_ids = take_i64(:mio_contact_pairs_master_cell, n)
        master_ids .= ifelse.(master_ids .< 0, -1, master_ids .+ 1)
        local_coordinates = reshape(take_f64(:mio_contact_pairs_local_coordinates, 3n), 3, n)
        closest_point = reshape(take_f64(:mio_contact_pairs_closest_point, 3n), 3, n)
        normal = reshape(take_f64(:mio_contact_pairs_normal, 3n), 3, n)
        (slave_point=slave_ids, master_cell=master_ids,
         master_facet=take_i64(:mio_contact_pairs_master_facet, n),
         master_subfacet=take_i64(:mio_contact_pairs_master_subfacet, n),
         local_coordinates=local_coordinates, closest_point=closest_point,
         gap=take_f64(:mio_contact_pairs_gap, n), normal=normal,
         unmatched=take_i64(:mio_contact_pairs_unmatched, nu) .+ 1)
    finally
        ccall(_sym(:mio_contact_pairs_result_free), Cvoid, (Ptr{Cvoid},), handle)
    end
end

"""
    split_interface(mesh, side; add_cohesive=false) -> (; mesh, num_duplicated_points, num_cohesive_cells)

Duplicate incident-cell fans across a Side region while preserving point/cell
data and regions. Polyhedra are rejected. See `doc/region_adjacency.md`.
"""
function split_interface(m::Mesh, side::AbstractString; add_cohesive::Bool=false)
    side_c = Vector{UInt8}(codeunits(String(side) * "\0"))
    return GC.@preserve side_c begin
        selector = _CRegionSelector(Cstring(pointer(side_c)), Int32(2), Int32(0), Int64(-2),
                                    Int64(-2), (Int64(0), Int64(0)))
        opts = _CSplitInterfaceOpts(Int32(add_cohesive), Int32(0), ntuple(_ -> Int64(0), 4))
        duplicated, cohesive = Ref{Int64}(0), Ref{Int64}(0)
        ptr = ccall(_sym(:mio_split_interface), Ptr{Cvoid},
                    (Ptr{Cvoid}, Ref{_CRegionSelector}, Ref{_CSplitInterfaceOpts},
                     Ptr{Int64}, Ptr{Int64}), _handle(m), Ref(selector), Ref(opts),
                    duplicated, cohesive)
        mesh = Mesh(_check_ptr(ptr))
        (mesh=mesh, num_duplicated_points=Int(duplicated[]),
         num_cohesive_cells=Int(cohesive[]))
    end
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
