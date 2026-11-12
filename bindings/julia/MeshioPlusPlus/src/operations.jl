# The operations layer: computations *on* a mesh, not file formats. Full parity
# with the Fortran module's `type(mio_mesh)` procedures.
#
# Ownership rule for the opaque C results (`mio_split_result`,
# `mio_partition_result`, `mio_reorder_result`, `mio_refine_result`,
# `mio_convert_cells_result`, `mio_decimate_result`): read the maps and
# counters first, then TAKE the mesh(es) out of the result, then free it. Every
# `Mesh` handed to the caller therefore owns its handle, and no Julia value can
# dangle when the result goes away.

"""
Copy an int64 index map, shifting to 1-based. The C API's `-1` sentinel
("pruned" / "absent" / "not in this piece") becomes `0`, which is never a valid
1-based index and so stays unambiguous. This is exactly the Fortran module's
rule; the two bindings must agree.
"""
function _shift_map(ptr::Ptr{Cvoid}, n::Integer)
    n <= 0 && return Int[]
    src = unsafe_wrap(Array, Ptr{Int64}(ptr), (Int(n),); own=false)
    out = Vector{Int}(undef, Int(n))
    @inbounds for i in eachindex(out)
        v = src[i]
        out[i] = v < 0 ? 0 : Int(v) + 1
    end
    out
end

"""Read a `(result) -> (data, dtype, n)` zero-copy map getter into a 1-based copy."""
function _result_map(result::Ptr{Cvoid}, sym::Symbol)
    data = Ref{Ptr{Cvoid}}(C_NULL); dt = Ref{Cint}(0); n = Ref{Int64}(0)
    _check(ccall(_sym(sym), Cint,
                 (Ptr{Cvoid}, Ptr{Ptr{Cvoid}}, Ptr{Cint}, Ptr{Int64}), result, data, dt, n))
    _shift_map(data[], n[])
end

"""Read a `(result, block) -> (data, dtype, n)` map getter into a 1-based copy."""
function _result_block_map(result::Ptr{Cvoid}, sym::Symbol, block::Integer)
    data = Ref{Ptr{Cvoid}}(C_NULL); dt = Ref{Cint}(0); n = Ref{Int64}(0)
    _check(ccall(_sym(sym), Cint,
                 (Ptr{Cvoid}, Int64, Ptr{Ptr{Cvoid}}, Ptr{Cint}, Ptr{Int64}),
                 result, Int64(block), data, dt, n))
    _shift_map(data[], n[])
end

"""Collect all per-block cell maps of an opaque result."""
function _result_cell_maps(result::Ptr{Cvoid}, count_sym::Symbol, map_sym::Symbol)
    n = _check_count(ccall(_sym(count_sym), Int64, (Ptr{Cvoid},), result), string(count_sym))
    [_result_block_map(result, map_sym, b) for b in 0:(n-1)]
end

# --- surface / skin / quality ------------------------------------------------

"""
    extract_surface(mesh; record_parent_ids=false) -> Mesh

Boundary of the mesh's highest-dimension cells: a volume mesh yields its
boundary faces, a 2-D surface mesh its boundary edges.

`record_parent_ids` attaches an Int64 `"surface:parent_cell"` cell-data array
naming each output facet's owning input cell.

Named regions are **not** carried through (the operation warns); see
`doc/extract_surface.md`.
"""
extract_surface(m::Mesh; record_parent_ids::Bool=false) =
    Mesh(_check_ptr(ccall(_sym(:mio_extract_surface), Ptr{Cvoid}, (Ptr{Cvoid}, Cint),
                          _handle(m), record_parent_ids ? Cint(1) : Cint(0))))

"""
    extract_skin(mesh; linearize=false) -> Mesh

Boundary skin of a volume mesh as a new surface mesh. `linearize` emits corner
nodes only (triangle/quad output).
"""
extract_skin(m::Mesh; linearize::Bool=false) =
    Mesh(_check_ptr(ccall(_sym(:mio_extract_skin), Ptr{Cvoid}, (Ptr{Cvoid}, Cint),
                          _handle(m), linearize ? Cint(1) : Cint(0))))

"""
    attach_quality(mesh) -> Mesh

A copy of the mesh with the per-cell quality metrics attached as cell data
(names `"quality:<metric>"`). See `doc/mesh_quality.md`.
"""
attach_quality(m::Mesh) =
    Mesh(_check_ptr(ccall(_sym(:mio_attach_quality), Ptr{Cvoid}, (Ptr{Cvoid},), _handle(m))))

"""
    quality_counts(mesh) -> (; num_cells, num_inverted, num_degenerate)
"""
function quality_counts(m::Mesh)
    nc = Ref{Int64}(0); ni = Ref{Int64}(0); nd = Ref{Int64}(0)
    _check(ccall(_sym(:mio_quality_counts), Cint,
                 (Ptr{Cvoid}, Ptr{Int64}, Ptr{Int64}, Ptr{Int64}), _handle(m), nc, ni, nd))
    (num_cells=Int(nc[]), num_inverted=Int(ni[]), num_degenerate=Int(nd[]))
end

# --- geometry-preserving transforms -----------------------------------------

"""
    transform(mesh, matrix; rotate_vector_data=false) -> Mesh

Affine transform of the point coordinates. `matrix` is 4×4; a point `p` maps to
`M * [p.x, p.y, p.z, 1]`. Connectivity and data are carried through; only the
points change, unless `rotate_vector_data` also rotates vector (trailing
dimension 3) and tensor (trailing dimension 9) point data by the transform's
3×3 linear block.

The matrix is passed to the C API row-major, so a `4×4` Julia matrix is
transposed on the way out — this is the one place a transpose is correct,
because a 4×4 affine matrix is a mathematical object, not a mesh array.
"""
function transform(m::Mesh, matrix::AbstractMatrix{<:Real};
                   rotate_vector_data::Bool=false)
    size(matrix) == (4, 4) || throw(ArgumentError(
        "the affine matrix must be 4x4, got $(size(matrix))"))
    rowmajor = Float64[matrix[i, j] for i in 1:4 for j in 1:4]
    ptr = GC.@preserve rowmajor ccall(_sym(:mio_transform), Ptr{Cvoid},
                                      (Ptr{Cvoid}, Ptr{Cdouble}, Cint),
                                      _handle(m), pointer(rowmajor),
                                      rotate_vector_data ? Cint(1) : Cint(0))
    Mesh(_check_ptr(ptr))
end

"""
    clean(mesh; weld=false, atol=1e-12, remove_orphans=true, drop_degenerate=true,
          drop_duplicate_cells=true)
        -> (; mesh, points_welded, points_removed_orphan,
              cells_dropped_degenerate, cells_dropped_duplicate)

One-pass weld / prune / de-duplicate. The default set matches the CLI's
no-flag behaviour: no welding, the other three steps on.
"""
function clean(m::Mesh; weld::Bool=false, atol::Real=1e-12, remove_orphans::Bool=true,
               drop_degenerate::Bool=true, drop_duplicate_cells::Bool=true)
    w = Ref{Int64}(0); o = Ref{Int64}(0); dg = Ref{Int64}(0); du = Ref{Int64}(0)
    ptr = ccall(_sym(:mio_clean), Ptr{Cvoid},
                (Ptr{Cvoid}, Cint, Cdouble, Cint, Cint, Cint,
                 Ptr{Int64}, Ptr{Int64}, Ptr{Int64}, Ptr{Int64}),
                _handle(m), weld ? Cint(1) : Cint(0), Float64(atol),
                remove_orphans ? Cint(1) : Cint(0), drop_degenerate ? Cint(1) : Cint(0),
                drop_duplicate_cells ? Cint(1) : Cint(0), w, o, dg, du)
    (mesh=Mesh(_check_ptr(ptr)), points_welded=Int(w[]),
     points_removed_orphan=Int(o[]), cells_dropped_degenerate=Int(dg[]),
     cells_dropped_duplicate=Int(du[]))
end

"""
    smooth(mesh; method="taubin", iterations=10, lambda=-1.0, mu=-0.34,
           fix_boundary=true, preserve_features=true, feature_angle=30.0,
           guard_inversion=true)
        -> (; mesh, num_nodes_moved, max_displacement, num_skipped_inversion)

Relax the point coordinates toward their edge-neighbour centroids. Topology and
every data value pass through unchanged — only the points move.

`method` is `"taubin"` (default, shrink-free), `"laplacian"` (stronger per
pass, shrinks the mesh), or `"odt"` (optimal-Delaunay-triangulation
smoothing — tet-only, raises on any other block; moves each free interior
vertex toward the volume-weighted average of its incident tets'
circumcenters).

A **negative `lambda` means "this method's own default"**: 0.5 for
`"laplacian"`, 0.33 for `"taubin"`, 0.9 for `"odt"`. Taubin additionally
needs `mu < -lambda < 0`; `mu` is silently ignored by `"laplacian"` and
`"odt"`. See `doc/smooth.md`.

`frozen` is an optional vector of **1-based** point ids to pin outright,
unioned with the boundary/feature pins; shifted to 0-based across the C ABI,
as every index array in this binding. An out-of-range id raises by name.
"""
function smooth(m::Mesh; method::AbstractString="taubin", iterations::Integer=10,
                lambda::Real=-1.0, mu::Real=-0.34, fix_boundary::Bool=true,
                preserve_features::Bool=true, feature_angle::Real=30.0,
                guard_inversion::Bool=true, frozen=nothing)
    moved = Ref{Int64}(0); disp = Ref{Cdouble}(0.0); skipped = Ref{Int64}(0)
    method_c = Vector{UInt8}(codeunits(String(method)))
    push!(method_c, 0x00)
    ids = frozen === nothing ? Int64[] : Vector{Int64}(vec(collect(frozen))) .- Int64(1)
    ptr = GC.@preserve ids method_c begin
        opts = _CSmoothOpts(Cstring(pointer(method_c)), Int32(iterations), Int32(0),
                            Cdouble(lambda), Cdouble(mu),
                            fix_boundary ? Int32(1) : Int32(0),
                            preserve_features ? Int32(1) : Int32(0),
                            Cdouble(feature_angle),
                            guard_inversion ? Int32(1) : Int32(0), Int32(0),
                            isempty(ids) ? Ptr{Int64}(C_NULL) : pointer(ids),
                            Int64(length(ids)), (Int64(0), Int64(0), Int64(0), Int64(0)))
        ccall(_sym(:mio_smooth_ex), Ptr{Cvoid},
             (Ptr{Cvoid}, Ref{_CSmoothOpts}, Ptr{Int64}, Ptr{Cdouble}, Ptr{Int64}),
             _handle(m), Ref(opts), moved, disp, skipped)
    end
    (mesh=Mesh(_check_ptr(ptr)), num_nodes_moved=Int(moved[]),
     max_displacement=Float64(disp[]), num_skipped_inversion=Int(skipped[]))
end

"""
    optimize_volume(mesh; max_iterations=10, relocate=true, flip=true,
                    preserve_boundary=true, min_improvement=1e-6,
                    objective="worst", max_ring_size=7)
        -> (; mesh, num_flips, num_23_flips, num_32_flips, num_edge_removals,
              num_multi_face_removals, num_vertices_moved, num_worst_smoothed,
              num_tets, min_quality_before, min_quality_after,
              min_dihedral_before, min_dihedral_after)

ODT-remesh a tetrahedral mesh: raise its worst element quality by relocating
vertices AND changing connectivity, predicate-free. `objective="worst"` (the
default) raises the least dihedral angle and never lowers it -- guarded ODT
relocation, 2-3 flips, edge and multi-face removal over rings of up to
`max_ring_size` tets, worst-cell smoothing; `"energy"` is the ODT energy with
scaled-Jacobian flips. The genuine "ODT
remeshing" sibling of [`remesh_volume`](@ref) (which generates a fresh lattice
mesh) and of `smooth(method="odt")` (which only moves points on fixed
connectivity). Tet-only, raises on any other block. The point set is invariant,
so point data and named Point regions carry; cell data and Cell/Side regions are
dropped. With `preserve_boundary` the boundary surface is exactly preserved.
See `doc/optimize_volume.md`.
"""
function optimize_volume(m::Mesh; max_iterations::Integer=10, relocate::Bool=true,
                         flip::Bool=true, preserve_boundary::Bool=true,
                         min_improvement::Real=1e-6, objective::AbstractString="worst",
                         max_ring_size::Integer=7)
    obj = objective == "worst" ? Int32(0) : objective == "energy" ? Int32(1) :
          throw(ArgumentError("meshio++: optimize_volume: objective must be \"worst\" or \"energy\""))
    opts = _COptimizeVolumeOpts(Int32(max_iterations), Int32(relocate), Int32(flip),
                                Int32(preserve_boundary), Float64(min_improvement), obj,
                                Int32(max_ring_size), ntuple(_ -> Int64(0), 8))
    report = Ref{_COptimizeVolumeReport}()
    ptr = ccall(_sym(:mio_optimize_volume_ex), Ptr{Cvoid},
                (Ptr{Cvoid}, Ref{_COptimizeVolumeOpts}, Ptr{_COptimizeVolumeReport}),
                _handle(m), Ref(opts), report)
    r = _check_ptr(ptr)
    rep = report[]
    (mesh=Mesh(r), num_flips=Int(rep.num_flips), num_23_flips=Int(rep.num_23_flips),
     num_32_flips=Int(rep.num_32_flips), num_edge_removals=Int(rep.num_edge_removals),
     num_multi_face_removals=Int(rep.num_multi_face_removals),
     num_vertices_moved=Int(rep.num_vertices_moved),
     num_worst_smoothed=Int(rep.num_worst_smoothed), num_tets=Int(rep.num_tets),
     min_quality_before=rep.min_quality_before, min_quality_after=rep.min_quality_after,
     min_dihedral_before=rep.min_dihedral_before, min_dihedral_after=rep.min_dihedral_after)
end

const _ADAPT_VOLUME_QUALITIES = Dict("min_dihedral" => Int32(0), "mean_ratio" => Int32(1),
                                     "amips" => Int32(2))

# Fill a mio_adapt_volume_opts from keyword values (`nothing` keeps the default).
function _adapt_volume_opts!(opts::Ref{_CAdaptVolumeOpts}; hmin=nothing, hmax=nothing,
                             hgrad=nothing, hausd=nothing, hausd_relative=nothing,
                             curvature_sizing=nothing, feature_angle=nothing,
                             corner_angle=nothing, max_normal_deviation=nothing,
                             preserve_boundary=nothing, preserve_features=nothing,
                             quality_metric=nothing, quality_min=nothing,
                             max_iterations=nothing, target_min_dihedral=nothing,
                             collapse=nothing, split=nothing, swap=nothing, smooth=nothing,
                             boundary_sizing=nothing)
    ccall(_sym(:mio_adapt_volume_opts_init), Cvoid, (Ref{_CAdaptVolumeOpts},), opts)
    scalars = (; hmin, hmax, hgrad, hausd, feature_angle, corner_angle, max_normal_deviation,
               quality_min, target_min_dihedral, max_iterations)
    for (k, v) in pairs(scalars)
        v === nothing || _setopt!(opts, k, v)
    end
    flags = (; hausd_relative, curvature_sizing, preserve_boundary, preserve_features,
             collapse, split, swap, smooth)
    for (k, v) in pairs(flags)
        v === nothing || _setopt!(opts, k, Bool(v) ? Int32(1) : Int32(0))
    end
    if quality_metric !== nothing
        code = get(_ADAPT_VOLUME_QUALITIES, String(quality_metric), nothing)
        code === nothing && throw(ArgumentError(
            "meshio++: adapt_volume: quality_metric must be \"min_dihedral\", \"mean_ratio\" or \"amips\""))
        _setopt!(opts, :quality_metric, code)
    end
    if boundary_sizing !== nothing
        b = String(boundary_sizing)
        b == "input" || b == "geometry" || throw(ArgumentError(
            "meshio++: adapt_volume: boundary_sizing must be \"input\" or \"geometry\""))
        _setopt!(opts, :boundary_sizing, b == "geometry" ? Int32(1) : Int32(0))
    end
    opts
end

_adapt_volume_named(rep::_CAdaptVolumeReport) =
    (num_tets_before=Int(rep.num_tets_before), num_tets_after=Int(rep.num_tets_after),
     num_points_before=Int(rep.num_points_before), num_points_after=Int(rep.num_points_after),
     num_collapsed=Int(rep.num_collapsed), num_split=Int(rep.num_split),
     num_swapped=Int(rep.num_swapped), num_smoothed=Int(rep.num_smoothed),
     num_rejected_by_link=Int(rep.num_rejected_by_link),
     num_rejected_by_inversion=Int(rep.num_rejected_by_inversion),
     num_rejected_by_quality=Int(rep.num_rejected_by_quality),
     num_rejected_by_length=Int(rep.num_rejected_by_length),
     num_rejected_by_class=Int(rep.num_rejected_by_class),
     num_rejected_by_hausdorff=Int(rep.num_rejected_by_hausdorff),
     num_rejected_by_normal=Int(rep.num_rejected_by_normal),
     num_corners=Int(rep.num_corners), num_ridge_vertices=Int(rep.num_ridge_vertices),
     num_chains=Int(rep.num_chains), num_patches=Int(rep.num_patches),
     iterations=Int(rep.iterations), min_dihedral_before=rep.min_dihedral_before,
     min_dihedral_after=rep.min_dihedral_after, max_dihedral_before=rep.max_dihedral_before,
     max_dihedral_after=rep.max_dihedral_after, mean_quality_before=rep.mean_quality_before,
     mean_quality_after=rep.mean_quality_after,
     max_boundary_deviation=rep.max_boundary_deviation, hmin=rep.hmin, hmax=rep.hmax,
     hausd=rep.hausd, volume_before=rep.volume_before, volume_after=rep.volume_after)

"""
    adapt_volume(m; hmin, hmax, hgrad, hausd, hausd_relative, size_array,
                 curvature_sizing, feature_angle, corner_angle, max_normal_deviation,
                 preserve_boundary, preserve_features, quality_metric, quality_min,
                 max_iterations, target_min_dihedral, collapse, split, swap, smooth)

Remesh a tetrahedral mesh to a graded isotropic size field: collapse, split, swap
and smooth under a boundary Hausdorff bound, keeping corners, ridges, materials,
regions and boundary tags. Every keyword defaults to `mio_adapt_volume_opts_init`'s
value. Returns `(; mesh, num_tets_before, num_tets_after, ...)` with every field of
the C report (see doc/adapt_volume.md).
"""
function adapt_volume(m::Mesh; size_array=nothing, kwargs...)
    opts = Ref{_CAdaptVolumeOpts}()
    _adapt_volume_opts!(opts; kwargs...)
    name = size_array === nothing ? "" : String(size_array)
    report = Ref{_CAdaptVolumeReport}()
    ptr = GC.@preserve name begin
        isempty(name) || _setopt!(opts, :size_array, Base.unsafe_convert(Cstring, name))
        ccall(_sym(:mio_adapt_volume), Ptr{Cvoid},
              (Ptr{Cvoid}, Ref{_CAdaptVolumeOpts}, Ref{_CAdaptVolumeReport}),
              _handle(m), opts, report)
    end
    r = _check_ptr(ptr)
    Base.merge((mesh=Mesh(r),), _adapt_volume_named(report[]))
end

# --- subsetting --------------------------------------------------------------

_crop_mode(mode::Symbol) = mode === :all ? Cint(0) : mode === :any ? Cint(1) :
                           throw(ArgumentError("mode must be :all or :any, got :$mode"))

"""
    crop_bbox(mesh, lo, hi; mode=:all, record_ids=false) -> Mesh

Keep the cells inside the axis-aligned box `lo .<= p .<= hi` (3-element
vectors). `mode` is `:all` (every node inside) or `:any` (at least one).
`record_ids` attaches `crop:original_point_id` / `crop:original_cell_id`.
"""
function crop_bbox(m::Mesh, lo::AbstractVector{<:Real}, hi::AbstractVector{<:Real};
                   mode::Symbol=:all, record_ids::Bool=false)
    l = Float64[lo[1], lo[2], lo[3]]; h = Float64[hi[1], hi[2], hi[3]]
    ptr = GC.@preserve l h ccall(_sym(:mio_crop_bbox), Ptr{Cvoid},
                                 (Ptr{Cvoid}, Ptr{Cdouble}, Ptr{Cdouble}, Cint, Cint),
                                 _handle(m), pointer(l), pointer(h), _crop_mode(mode),
                                 record_ids ? Cint(1) : Cint(0))
    Mesh(_check_ptr(ptr))
end

"""
    crop_plane(mesh, point, normal; mode=:all, record_ids=false) -> Mesh

Keep the cells in the half-space `(p - point) . normal >= 0`. Unlike
[`slice`](@ref), whole cells are kept and the dimension is unchanged.
"""
function crop_plane(m::Mesh, point::AbstractVector{<:Real},
                    normal::AbstractVector{<:Real}; mode::Symbol=:all,
                    record_ids::Bool=false)
    p = Float64[point[1], point[2], point[3]]; n = Float64[normal[1], normal[2], normal[3]]
    ptr = GC.@preserve p n ccall(_sym(:mio_crop_plane), Ptr{Cvoid},
                                 (Ptr{Cvoid}, Ptr{Cdouble}, Ptr{Cdouble}, Cint, Cint),
                                 _handle(m), pointer(p), pointer(n), _crop_mode(mode),
                                 record_ids ? Cint(1) : Cint(0))
    Mesh(_check_ptr(ptr))
end

"""
    crop_predicate(mesh, array; compare="<", value=0.0, record_ids=false) -> Mesh

Crop to the cells whose value in the scalar `cell_data` array `array` satisfies
a comparison -- one of `"<"`, `"<="`, `">"`, `">="`, `"=="`, `"!="`, the same
vocabulary [`refine`](@ref)'s `where_op` uses, evaluated by the same C++
function. **A non-finite cell value never matches**, whatever the comparison.

There is deliberately no `mode`: `crop_bbox`/`crop_plane` test points and then
need an all/any rule, whereas a `cell_data` predicate is already one value per
cell and has nothing to reduce.

Inside/outside a surface composes rather than being its own function:

```julia
f = distance_to_surface(m, skin; location=:center).mesh
inside = crop_predicate(f, "sdf:distance"; compare="<", value=0.0)
```
"""
function crop_predicate(m::Mesh, array::AbstractString; compare="<", value::Real=0.0,
                        record_ids::Bool=false)
    code = get(_REFINE_COMPARES, String(compare)) do
        throw(ArgumentError("meshio++: unknown comparison '$(compare)'"))
    end
    ptr = ccall(_sym(:mio_crop_predicate), Ptr{Cvoid},
                (Ptr{Cvoid}, Cstring, Cint, Cdouble, Cint),
                _handle(m), String(array), code, Cdouble(value),
                record_ids ? Cint(1) : Cint(0))
    Mesh(_check_ptr(ptr))
end

"""
    slice(mesh, origin, normal; record_parent_ids=false) -> Mesh

Planar cross-section: the actual **intersection** with the plane, one
topological dimension below the cut cells (a volume mesh yields a
triangle/quad surface, a 2-D surface mesh a line mesh). Contrast
[`crop_plane`](@ref), which keeps whole cells on one side.

Empty when the plane misses the geometry — not an error. See `doc/slice.md`.
"""
function slice(m::Mesh, origin::AbstractVector{<:Real}, normal::AbstractVector{<:Real};
               record_parent_ids::Bool=false, chain::Bool=false, min_length::Real=0.0)
    o = Float64[origin[1], origin[2], origin[3]]
    n = Float64[normal[1], normal[2], normal[3]]
    ptr = GC.@preserve o n ccall(_sym(:mio_slice_ex), Ptr{Cvoid},
                                 (Ptr{Cvoid}, Ptr{Cdouble}, Ptr{Cdouble}, Cint, Cint, Cdouble),
                                 _handle(m), pointer(o), pointer(n),
                                 record_parent_ids ? Cint(1) : Cint(0),
                                 chain ? Cint(1) : Cint(0), Cdouble(min_length))
    Mesh(_check_ptr(ptr))
end

"""
    clip(m, origin, normal; keep_below=true, cap=true)
        -> (; mesh, num_kept, num_split, num_cap_triangles, num_loops, num_open_chains, warnings)

Cut a surface by the plane through `origin` normal to `normal`, keep a side and
cap the cut's closed loops (roadmap 7.8.3). See `doc/clip.md`.
"""
function clip(m::Mesh, origin, normal; keep_below::Bool=true, cap::Bool=true)
    o = Cdouble[origin[1], origin[2], origin[3]]
    n = Cdouble[normal[1], normal[2], normal[3]]
    rep = Ref{_CClipReport}()
    ptr = GC.@preserve o n ccall(_sym(:mio_clip), Ptr{Cvoid},
                                 (Ptr{Cvoid}, Ptr{Cdouble}, Ptr{Cdouble}, Int32, Int32,
                                  Ref{_CClipReport}), _handle(m), o, n, Int32(keep_below),
                                 Int32(cap), rep)
    out = Mesh(_check_ptr(ptr))
    r = rep[]
    text = _getstring((buf, k) -> ccall(_sym(:mio_clip_warnings), Int64, (Ptr{UInt8}, Int64), buf, k))
    (mesh=out, num_kept=Int(r.num_kept), num_split=Int(r.num_split),
     num_cap_triangles=Int(r.num_cap_triangles), num_loops=Int(r.num_loops),
     num_open_chains=Int(r.num_open_chains),
     warnings=isempty(text) ? String[] : String.(Base.split(text, '\n')))
end

"""
    profile_solid(sketch, features; step_path="", scan=nothing, provider="", deflection=0.0)
        -> (; mesh, volume, area, num_faces, valid, warnings)

Build a `fit_curves` sketch into a solid by `features` (text,
`;`-separated: `extrude:D[:symmetric]`, `revolve:OX,OY,OZ:DX,DY,DZ[:ANGLE]`,
`fillet:R`, `chamfer:D`, `shell:W`) through a CAD provider (OpenCASCADE;
roadmap 7.7.4). See `doc/profile_solid.md`.
"""
function profile_solid(sketch::Mesh, features::AbstractString; step_path::AbstractString="",
                       scan::Union{Mesh,Nothing}=nothing, provider::AbstractString="",
                       deflection::Real=0.0)
    text = String(features)
    n = ccall(_sym(:mio_profile_features_from_text), Int64,
              (Cstring, Ptr{_CProfileFeature}, Int64), text, C_NULL, 0)
    n < 0 && _check_ptr(C_NULL)
    feats = Vector{_CProfileFeature}(undef, max(n, 1))
    ccall(_sym(:mio_profile_features_from_text), Int64,
          (Cstring, Ptr{_CProfileFeature}, Int64), text, feats, n)
    rep = Ref{_CProfileSolidReport}()
    ptr = ccall(_sym(:mio_profile_solid), Ptr{Cvoid},
                (Ptr{Cvoid}, Ptr{_CProfileFeature}, Int64, Cstring, Cstring, Cdouble, Ptr{Cvoid},
                 Ref{_CProfileSolidReport}),
                _handle(sketch), feats, n, String(step_path), String(provider),
                Float64(deflection), scan === nothing ? C_NULL : _handle(scan), rep)
    out = Mesh(_check_ptr(ptr))
    r = rep[]
    text_w = _getstring((buf, k) -> ccall(_sym(:mio_profile_solid_warnings), Int64,
                                          (Ptr{UInt8}, Int64), buf, k))
    (mesh=out, volume=r.volume, area=r.area, num_faces=Int(r.num_faces), valid=r.valid != 0,
     warnings=isempty(text_w) ? String[] : String.(Base.split(text_w, '\n')))
end

const _RELATIONS = ("coaxial", "parallel", "perpendicular", "equal_radius")

# (a, b, relation name) tuples as int64 triples.
function _relation_triples(list)
    out = Int64[]
    for (a, b, kind) in list
        k = findfirst(==(String(kind)), _RELATIONS)
        k === nothing && throw(ArgumentError(
            "meshio++: constrain_primitives: a relation is coaxial, parallel, perpendicular or equal_radius"))
        append!(out, (Int64(a), Int64(b), Int64(k - 1)))
    end
    out
end

"""
    constrain_primitives(m; face_tolerance=0.0, angle_tolerance=-1.0, distance_tolerance=-1.0,
                         min_fraction=-1.0, force=(), forbid=())
        -> (; mesh, num_constrained, num_groups, face_tolerance, relations, unmet)

Fit related primitives of a `fit_primitives` skin together under control
(roadmap 7.7.3): `force` and `forbid` hold `(patch_a, patch_b, relation)`
tuples, patches by `segment:patch` label. `relations` are named tuples
`(; patch_a, patch_b, relation, angle_deg, distance, forced)`. See
`doc/constrain_primitives.md`.
"""
function constrain_primitives(m::Mesh; face_tolerance::Real=0.0, angle_tolerance::Real=-1.0,
                              distance_tolerance::Real=-1.0, min_fraction::Real=-1.0,
                              force=(), forbid=())
    fo = _relation_triples(force)
    fb = _relation_triples(forbid)
    rep = Ref{_CConstrainPrimitivesReport}()
    ptr = GC.@preserve fo fb begin
        opts = _CConstrainPrimitivesOpts(Cdouble(face_tolerance), Cdouble(angle_tolerance),
                                         Cdouble(distance_tolerance), Cdouble(min_fraction),
                                         isempty(fo) ? Ptr{Int64}(C_NULL) : pointer(fo),
                                         Int64(length(fo) ÷ 3),
                                         isempty(fb) ? Ptr{Int64}(C_NULL) : pointer(fb),
                                         Int64(length(fb) ÷ 3), (0, 0, 0, 0))
        ccall(_sym(:mio_constrain_primitives), Ptr{Cvoid},
              (Ptr{Cvoid}, Ref{_CConstrainPrimitivesOpts}, Ref{_CConstrainPrimitivesReport}),
              _handle(m), Ref(opts), rep)
    end
    out = Mesh(_check_ptr(ptr))
    r = rep[]
    relations = map(0:r.num_relations-1) do k
        a = Ref{Int64}(0); b = Ref{Int64}(0); kind = Ref{Int32}(0)
        ang = Ref{Cdouble}(0); dst = Ref{Cdouble}(0); fc = Ref{Int32}(0)
        ccall(_sym(:mio_constrain_primitives_relation), Cint,
              (Int64, Ref{Int64}, Ref{Int64}, Ref{Int32}, Ref{Cdouble}, Ref{Cdouble}, Ref{Int32}),
              k, a, b, kind, ang, dst, fc)
        (; patch_a=Int(a[]), patch_b=Int(b[]), relation=_RELATIONS[kind[] + 1],
         angle_deg=ang[], distance=dst[], forced=fc[] != 0)
    end
    unmet = map(0:r.num_unmet-1) do k
        a = Ref{Int64}(0); b = Ref{Int64}(0); kind = Ref{Int32}(0)
        ccall(_sym(:mio_constrain_primitives_unmet), Cint,
              (Int64, Ref{Int64}, Ref{Int64}, Ref{Int32}), k, a, b, kind)
        (Int(a[]), Int(b[]), _RELATIONS[kind[] + 1])
    end
    (mesh=out, num_constrained=Int(r.num_constrained), num_groups=Int(r.num_groups),
     face_tolerance=r.face_tolerance, relations=relations, unmet=unmet)
end

"""
    face_support(solid, scan; tolerance=0.0, measured_share=0.9, inferred_share=0.1,
                 hole_data="repair:hole")
        -> (; mesh, num_measured, num_partial, num_inferred, num_faces, share_within, tolerance)

Class each face of a `fit_cad` re-mesh by how much of it the scan supports
(roadmap 7.7.1). See `doc/fit_cad_loop.md`.
"""
function face_support(solid::Mesh, scan::Mesh; tolerance::Real=0.0, measured_share::Real=0.9,
                      inferred_share::Real=0.1, hole_data::AbstractString="repair:hole")
    rep = Ref{_CFaceSupportReport}()
    hd = String(hole_data)
    ptr = GC.@preserve hd begin
        opts = _CFaceSupportOpts(Cdouble(tolerance), Cdouble(measured_share),
                                 Cdouble(inferred_share), Base.unsafe_convert(Cstring, hd),
                                 (0, 0, 0, 0))
        ccall(_sym(:mio_face_support), Ptr{Cvoid},
              (Ptr{Cvoid}, Ptr{Cvoid}, Ref{_CFaceSupportOpts}, Ref{_CFaceSupportReport}),
              _handle(solid), _handle(scan), Ref(opts), rep)
    end
    out = Mesh(_check_ptr(ptr))
    r = rep[]
    (mesh=out, num_measured=Int(r.num_measured), num_partial=Int(r.num_partial),
     num_inferred=Int(r.num_inferred), num_faces=Int(r.num_faces),
     share_within=r.share_within, tolerance=r.tolerance)
end

const _PART_CLASSES = ("any", "analytic", "cast", "freeform")
const _STOP_REASONS = ("target", "no improvement", "rounds", "nothing to refine")

"""
    fit_cad_iterate(scan; step_path="", part_class="any", target_within=0.0,
                    target_explained=0.0, segment_tolerance=0.0, fit_tolerance=0.0,
                    max_rounds=4, worst_faces=3, refine_factor=0.25, backend="auto")
        -> (; mesh, num_rounds, best_round, target_met, stop_reason, target_within,
              share_within, explained, hausdorff_max, num_faces, rounds)

Segment, fit and build a solid in rounds towards a target (roadmap 7.7.2);
`rounds` are named tuples `(; share_within, explained, hausdorff_max,
num_faces, num_resegmented)`. See `doc/fit_cad_loop.md`.
"""
function fit_cad_iterate(scan::Mesh; step_path::AbstractString="", part_class::AbstractString="any",
                         target_within::Real=0.0, target_explained::Real=0.0,
                         segment_tolerance::Real=0.0, fit_tolerance::Real=0.0,
                         max_rounds::Integer=4, worst_faces::Integer=3, refine_factor::Real=0.25,
                         backend::AbstractString="auto")
    pc = findfirst(==(String(part_class)), _PART_CLASSES)
    pc === nothing && throw(ArgumentError(
        "meshio++: fit_cad_iterate: part_class must be any, analytic, cast or freeform"))
    be = findfirst(==(String(backend)), ("auto", "occt", "native"))
    be === nothing && throw(ArgumentError("meshio++: fit_cad_iterate: backend must be auto, occt or native"))
    opts = _CFitCadIterateOpts(Cdouble(segment_tolerance), Cdouble(fit_tolerance),
                               Cdouble(target_within), Cdouble(target_explained),
                               Cdouble(refine_factor), Int32(pc - 1), Int32(be - 1),
                               Int64(max_rounds), Int64(worst_faces), Int64(1), (0, 0, 0, 0, 0))
    rep = Ref{_CFitCadIterateReport}()
    ptr = ccall(_sym(:mio_fit_cad_iterate), Ptr{Cvoid},
                (Ptr{Cvoid}, Cstring, Ref{_CFitCadIterateOpts}, Ref{_CFitCadIterateReport}),
                _handle(scan), String(step_path), Ref(opts), rep)
    out = Mesh(_check_ptr(ptr))
    r = rep[]
    rounds = map(0:r.num_rounds-1) do k
        sw = Ref{Cdouble}(0); ex = Ref{Cdouble}(0); hm = Ref{Cdouble}(0)
        nf = Ref{Int64}(0); nr = Ref{Int64}(0)
        ccall(_sym(:mio_fit_cad_iterate_round), Cint,
              (Int64, Ref{Cdouble}, Ref{Cdouble}, Ref{Cdouble}, Ref{Int64}, Ref{Int64}),
              k, sw, ex, hm, nf, nr)
        (; share_within=sw[], explained=ex[], hausdorff_max=hm[], num_faces=Int(nf[]),
         num_resegmented=Int(nr[]))
    end
    (mesh=out, num_rounds=Int(r.num_rounds), best_round=Int(r.best_round),
     target_met=r.target_met != 0, stop_reason=_STOP_REASONS[r.stop_reason + 1],
     target_within=r.target_within, share_within=r.share_within, explained=r.explained,
     hausdorff_max=r.hausdorff_max, num_faces=Int(r.num_faces), rounds=rounds)
end

# The shared tail of offset_surface, thicken and hollow.
function _offset_result(ptr, rep)
    out = Mesh(_check_ptr(ptr))
    r = rep[]
    text = _getstring((buf, k) -> ccall(_sym(:mio_offset_warnings), Int64, (Ptr{UInt8}, Int64), buf, k))
    (mesh=out, num_not_visible=Int(r.num_not_visible), num_clamped=Int(r.num_clamped),
     num_folded=Int(r.num_folded), num_rim_quads=Int(r.num_rim_quads), max_move=r.max_move,
     warnings=isempty(text) ? String[] : String.(Base.split(text, '\n')))
end

"""
    offset_surface(m, distance; smooth_iterations=0)
        -> (; mesh, num_not_visible, num_clamped, num_folded, num_rim_quads, max_move, warnings)

Move each point of a triangle and quad surface along its offset direction so
every face's plane moves by `distance` (roadmap 7.8.4). See `doc/offset_surface.md`.
"""
function offset_surface(m::Mesh, distance::Real; smooth_iterations::Integer=0)
    rep = Ref{_COffsetReport}()
    ptr = ccall(_sym(:mio_offset_surface), Ptr{Cvoid}, (Ptr{Cvoid}, Cdouble, Int32, Ref{_COffsetReport}),
                _handle(m), Float64(distance), Int32(smooth_iterations), rep)
    _offset_result(ptr, rep)
end

"""
    thicken(m, thickness; side="inward", smooth_iterations=0) -> (; mesh, ..., num_rim_quads, ...)

A closed shell `thickness` thick from a surface; `side` is `"inward"`,
`"outward"` or `"both"` (roadmap 7.8.4). See `doc/offset_surface.md`.
"""
function thicken(m::Mesh, thickness::Real; side::AbstractString="inward",
                 smooth_iterations::Integer=0)
    s = side == "inward" ? Int32(0) : side == "outward" ? Int32(1) : side == "both" ? Int32(2) :
        throw(ArgumentError("meshio++: thicken: side must be inward, outward or both"))
    rep = Ref{_COffsetReport}()
    ptr = ccall(_sym(:mio_thicken), Ptr{Cvoid}, (Ptr{Cvoid}, Cdouble, Int32, Int32, Ref{_COffsetReport}),
                _handle(m), Float64(thickness), s, Int32(smooth_iterations), rep)
    _offset_result(ptr, rep)
end

"""
    hollow(m, wall; smooth_iterations=0) -> (; mesh, ...)

A closed surface with an inner, reversed wall `wall` inside it (roadmap 7.8.4).
See `doc/offset_surface.md`.
"""
function hollow(m::Mesh, wall::Real; smooth_iterations::Integer=0)
    rep = Ref{_COffsetReport}()
    ptr = ccall(_sym(:mio_hollow), Ptr{Cvoid}, (Ptr{Cvoid}, Cdouble, Int32, Ref{_COffsetReport}),
                _handle(m), Float64(wall), Int32(smooth_iterations), rep)
    _offset_result(ptr, rep)
end

"""
    flip_normals(m; region="") -> (; mesh, num_flipped)

Reverse the winding of every cell, or of the cell region `region` (roadmap
7.8.2), with no point moved. See `doc/flip_normals.md`.
"""
function flip_normals(m::Mesh; region::AbstractString="")
    n = Ref{Int64}(0)
    ptr = ccall(_sym(:mio_flip_normals), Ptr{Cvoid}, (Ptr{Cvoid}, Cstring, Ref{Int64}),
                _handle(m), String(region), n)
    (mesh=Mesh(_check_ptr(ptr)), num_flipped=Int(n[]))
end

"""
    remove_islands(m; connectivity="edge", min_diagonal=0.01, min_area=0.0, min_cells=0,
                   relative=true, dry_run=false)
        -> (; mesh, num_components, num_removed, num_cells_removed, num_points_removed)

Remove the small connected components of `m` (roadmap 7.8.1). See
`doc/remove_islands.md`.
"""
function remove_islands(m::Mesh; connectivity::AbstractString="edge", min_diagonal::Real=0.01,
                        min_area::Real=0.0, min_cells::Integer=0, relative::Bool=true,
                        dry_run::Bool=false)
    conn = connectivity == "edge" ? Int32(0) :
           connectivity == "point" ? Int32(1) :
           throw(ArgumentError("meshio++: remove_islands: connectivity must be edge or point"))
    opts = _CRemoveIslandsOpts(conn, Int32(relative), Cdouble(min_diagonal), Cdouble(min_area),
                               Int64(min_cells), Int32(dry_run), Int32(0), ntuple(_ -> Int64(0), 3))
    rep = Ref{_CRemoveIslandsReport}()
    ptr = ccall(_sym(:mio_remove_islands), Ptr{Cvoid},
                (Ptr{Cvoid}, Ref{_CRemoveIslandsOpts}, Ref{_CRemoveIslandsReport}), _handle(m),
                Ref(opts), rep)
    out = Mesh(_check_ptr(ptr))
    r = rep[]
    (mesh=out, num_components=Int(r.num_components), num_removed=Int(r.num_removed),
     num_cells_removed=Int(r.num_cells_removed), num_points_removed=Int(r.num_points_removed))
end

"""
    measure(m, spec::AbstractString)
        -> (; mesh, num_elements, num_dimensions, num_tolerances, num_fail, report_json, warnings)

Measure `m` by the JSON spec `spec` (roadmap 7.2): elements fitted from regions,
seeds or patches, constructed or literal; dimensions and GD&T tolerances. The
mesh comes back with its element table; `report_json` is the full report. See
`doc/measure.md`.
"""
function measure(m::Mesh, spec::AbstractString)
    rep = Ref{_CMeasureReport}()
    ptr = ccall(_sym(:mio_measure), Ptr{Cvoid}, (Ptr{Cvoid}, Cstring, Ref{_CMeasureReport}),
                _handle(m), String(spec), rep)
    out = Mesh(_check_ptr(ptr))
    r = rep[]
    report = _getstring((buf, n) -> ccall(_sym(:mio_measure_report_json), Int64,
                                          (Ptr{UInt8}, Int64), buf, n))
    text = _getstring((buf, n) -> ccall(_sym(:mio_measure_warnings), Int64,
                                        (Ptr{UInt8}, Int64), buf, n))
    (mesh=out, num_elements=Int(r.num_elements), num_dimensions=Int(r.num_dimensions),
     num_tolerances=Int(r.num_tolerances), num_fail=Int(r.num_fail), report_json=report,
     warnings=isempty(text) ? String[] : String.(Base.split(text, '\n')))
end

"""
    align_datums(m, primary, secondary; origin="") -> (; mesh, matrix)

The 3-2-1 datum alignment of a measured mesh (roadmap 7.3.4): `primary` sets
Z, `secondary` X, `origin` (a point-like element or a third plane) the origin.
`matrix` is the 4 x 4 `transform` applies.
"""
function align_datums(m::Mesh, primary::AbstractString, secondary::AbstractString;
                      origin::AbstractString="")
    mat = zeros(Cdouble, 16)
    ptr = ccall(_sym(:mio_align_datums), Ptr{Cvoid},
                (Ptr{Cvoid}, Cstring, Cstring, Cstring, Ptr{Cdouble}), _handle(m),
                String(primary), String(secondary), String(origin), mat)
    out = Mesh(_check_ptr(ptr))
    (mesh=out, matrix=[mat[4 * (a - 1) + b] for a in 1:4, b in 1:4])
end

"""
    measure_element_primitive(m, name) -> (; type, params)

The element `name` of a measured mesh as a fitted primitive (1 plane, 2 sphere,
3 cylinder, 4 cone) and its ten parameters, a deviation reference (roadmap 7.4.2).
"""
function measure_element_primitive(m::Mesh, name::AbstractString)
    t = Ref{Int32}(0)
    p = zeros(Cdouble, 10)
    _check(ccall(_sym(:mio_measure_element_primitive), Cint,
                 (Ptr{Cvoid}, Cstring, Ref{Int32}, Ptr{Cdouble}), _handle(m), String(name), t, p))
    (type=Int(t[]), params=p)
end

"""
    write_step_elements(path, m; patches=Int64[], segments=64)
        -> (; num_planes, num_spheres, num_cylinders, num_cones, num_tori, num_skipped, warnings)

Write the fitted primitives of `m` (`segment`, then `fit_primitives`) to a STEP
file as exact geometry: bounded planes, balls, lidded cylinders, frusta and tori
(roadmap 7.6.5). `patches` selects patch labels (all fitted ones when empty).
"""
function write_step_elements(path::AbstractString, m::Mesh; patches=Int64[], segments::Integer=64)
    ids = Int64[Int64(p) for p in patches]
    rep = Ref{_CStepElementsReport}()
    GC.@preserve ids begin
        _check(ccall(_sym(:mio_write_step_elements), Cint,
                     (Cstring, Ptr{Cvoid}, Ptr{Int64}, Int64, Int64, Ref{_CStepElementsReport}),
                     String(path), _handle(m), isempty(ids) ? Ptr{Int64}(C_NULL) : pointer(ids),
                     Int64(length(ids)), Int64(segments), rep))
    end
    r = rep[]
    text = _getstring((buf, n) -> ccall(_sym(:mio_write_step_elements_warnings), Int64,
                                        (Ptr{UInt8}, Int64), buf, n))
    (num_planes=Int(r.num_planes), num_spheres=Int(r.num_spheres),
     num_cylinders=Int(r.num_cylinders), num_cones=Int(r.num_cones), num_tori=Int(r.num_tori),
     num_skipped=Int(r.num_skipped),
     warnings=isempty(text) ? String[] : String.(Base.split(text, '\n')))
end

"""
    silhouette(m, normal=(0, 0, 1); clearance=0.0, round_angle=180.0, arc_segments=64,
               extrude=false, region="", merge_tolerance=1e-9)
        -> (; mesh, num_loops, num_outer, num_holes, num_triangles, area, perimeter,
              min_height, max_height, warnings)

The outline of a surface (or a volume mesh's skin) seen along `normal`, with its
holes: a chained line mesh `fit_curves` reads, offset by `clearance` (> 0 grows it)
or, with `extrude=true`, the closed surface of the prism over it (roadmap 7.6.3).
See `doc/silhouette.md`.
"""
function silhouette(m::Mesh, normal=(0.0, 0.0, 1.0); clearance::Real=0.0, round_angle::Real=180.0,
                    arc_segments::Integer=64, extrude::Bool=false, region::AbstractString="",
                    merge_tolerance::Real=1e-9)
    n = (Float64(normal[1]), Float64(normal[2]), Float64(normal[3]))
    region_c = Vector{UInt8}(codeunits(String(region) * "\0"))
    rep = Ref{_CSilhouetteReport}()
    ptr = GC.@preserve region_c begin
        opts = _CSilhouetteOpts(n, Cdouble(clearance), Cdouble(round_angle), Int64(arc_segments),
                                Int32(extrude), Int32(0), Cdouble(merge_tolerance),
                                Cstring(pointer(region_c)), ntuple(_ -> Int64(0), 4))
        ccall(_sym(:mio_silhouette), Ptr{Cvoid},
              (Ptr{Cvoid}, Ref{_CSilhouetteOpts}, Ref{_CSilhouetteReport}), _handle(m), Ref(opts),
              rep)
    end
    out = Mesh(_check_ptr(ptr))
    r = rep[]
    text = _getstring((buf, n) -> ccall(_sym(:mio_silhouette_warnings), Int64,
                                        (Ptr{UInt8}, Int64), buf, n))
    (mesh=out, num_loops=Int(r.num_loops), num_outer=Int(r.num_outer),
     num_holes=Int(r.num_holes), num_triangles=Int(r.num_triangles), area=r.area,
     perimeter=r.perimeter, min_height=r.min_height, max_height=r.max_height,
     warnings=isempty(text) ? String[] : String.(Base.split(text, '\n')))
end

"""
    fit_curves(sections; tolerance=0.0, snap_angle=0.0, normal=nothing, min_arc_points=4,
               min_arc_sweep=5.0)
        -> (; mesh, pieces, frames, num_lines, num_arcs, num_circles, max_residual, tolerance,
              warnings)

Lines, arcs and circles from the chains of a chained line mesh (roadmap 7.6.2):
each piece is `(; type, start, stop, centre, radius, sweep, residual, num_points,
chain, frame)` in its frame's 2-D coordinates (`chain` and `frame` 0-based);
`frames` are 4 x 4 (columns e1, e2, normal, origin). See `doc/fit_curves.md`.
"""
function fit_curves(sections::Mesh; tolerance::Real=0.0, snap_angle::Real=0.0, normal=nothing,
                    min_arc_points::Integer=4, min_arc_sweep::Real=5.0)
    n = normal === nothing ? (0.0, 0.0, 0.0) : (Float64(normal[1]), Float64(normal[2]), Float64(normal[3]))
    opts = _CFitCurvesOpts(Cdouble(tolerance), Cdouble(snap_angle), n, Int64(min_arc_points),
                           Cdouble(min_arc_sweep), ntuple(_ -> Int64(0), 4))
    rep = Ref{_CFitCurvesReport}()
    ptr = ccall(_sym(:mio_fit_curves), Ptr{Cvoid},
                (Ptr{Cvoid}, Ref{_CFitCurvesOpts}, Ref{_CFitCurvesReport}), _handle(sections),
                Ref(opts), rep)
    out = Mesh(_check_ptr(ptr))
    r = rep[]
    names = ("line", "arc", "circle")
    pieces = map(0:(r.num_pieces - 1)) do i
        p = Ref{_CSketchPiece}()
        _check(ccall(_sym(:mio_fit_curves_piece), Cint, (Int64, Ref{_CSketchPiece}), i, p))
        q = p[]
        (type=names[q.type], start=collect(q.start), stop=collect(q.stop), centre=collect(q.centre),
         radius=q.radius, sweep=q.sweep, residual=q.residual, num_points=Int(q.num_points),
         chain=Int(q.chain), frame=Int(q.frame))
    end
    frames = map(0:(r.num_frames - 1)) do i
        f = zeros(Cdouble, 16)
        _check(ccall(_sym(:mio_fit_curves_frame), Cint, (Int64, Ptr{Cdouble}), i, f))
        [f[4 * (a - 1) + b] for a in 1:4, b in 1:4]
    end
    text = _getstring((buf, n) -> ccall(_sym(:mio_fit_curves_warnings), Int64,
                                        (Ptr{UInt8}, Int64), buf, n))
    (mesh=out, pieces=pieces, frames=frames, num_lines=Int(r.num_lines), num_arcs=Int(r.num_arcs),
     num_circles=Int(r.num_circles), max_residual=r.max_residual, tolerance=r.tolerance,
     warnings=isempty(text) ? String[] : String.(Base.split(text, '\n')))
end

"""
    sections(mesh, direction, offsets; min_length=0.0, record_parent_ids=false) -> Mesh

Chained sections of a surface at the parallel planes `n . x = offsets[k]` (`n`
the unit `direction`), merged into one line mesh with the `section:offset`,
`section:chain` and `section:closed` cell data (roadmap 7.6.1; chains 0-based).
See `doc/slice.md`.
"""
function sections(m::Mesh, direction::AbstractVector{<:Real}, offsets::AbstractVector{<:Real};
                  min_length::Real=0.0, record_parent_ids::Bool=false)
    d = Float64[direction[1], direction[2], direction[3]]
    offs = Float64.(collect(offsets))
    ptr = GC.@preserve d offs ccall(_sym(:mio_sections), Ptr{Cvoid},
                                    (Ptr{Cvoid}, Ptr{Cdouble}, Ptr{Cdouble}, Int64, Cdouble, Cint),
                                    _handle(m), pointer(d), pointer(offs), length(offs),
                                    Cdouble(min_length), record_parent_ids ? Cint(1) : Cint(0))
    Mesh(_check_ptr(ptr))
end

"""
    isosurface(mesh, array, isovalues; component=-1, record_parent_ids=false) -> Mesh

Level set of a scalar **point-data** field — the data-driven sibling of
[`slice`](@ref), sharing its marching-tetrahedra cutter.

`isovalues` may be a single number or a vector; all contours land in the one
returned mesh, tagged per cell with a Float64 `"iso:value"` and an Int64
`"iso:index"`. `component` selects a component of a multi-component array;
negative means the row magnitude.

Naming a *cell-data* array fails by name: cell data is piecewise constant and
has no level set — convert it with [`data_cell_to_point`](@ref) first. See
`doc/isosurface.md`.
"""
function isosurface(m::Mesh, array::AbstractString, isovalues;
                    component::Integer=-1, record_parent_ids::Bool=false)
    vals = Float64[Float64(v) for v in (isovalues isa Number ? (isovalues,) : isovalues)]
    isempty(vals) && throw(ArgumentError("at least one isovalue is required"))
    ptr = GC.@preserve vals ccall(_sym(:mio_isosurface), Ptr{Cvoid},
                                  (Ptr{Cvoid}, Cstring, Ptr{Cdouble}, Cint, Cint, Cint),
                                  _handle(m), array, pointer(vals), Cint(length(vals)),
                                  Cint(component), record_parent_ids ? Cint(1) : Cint(0))
    Mesh(_check_ptr(ptr))
end

"""
    gradient(mesh, array; operator=:gradient, method=:green_gauss, location=:cell,
             output="", component=-1, overwrite=false)
        -> (; mesh, num_skipped, num_fallback)

Differentiate a **point-data** field: its gradient, divergence or curl.

`operator` is `:gradient`, `:divergence` or `:curl`; `method` is `:green_gauss`
(exact for a linear field on any cell) or `:least_squares` (a fit over the
node-sharing neighbours); `location` is `:cell` or `:point`. The result is named
`"<array>:<operator>"` unless `output` overrides it.

An `nc`-component input yields `3 * nc` gradient components, flat and row-major
as `[component i][derivative j]`; divergence yields 1 and curl always 3, both
needing a 2- or 3-component field. `component` selects a single component of the
gradient; **negative means every component** — deliberately the opposite of
[`isosurface`](@ref), where negative means the row magnitude.

Cells that cannot be differentiated yield `NaN` and are counted in
`num_skipped`; least-squares cells with a degenerate neighbourhood fall back to
Green-Gauss and are counted in `num_fallback`.

Naming a *cell-data* array fails by name: cell data is piecewise constant and
has no derivative — convert it with [`data_cell_to_point`](@ref) first. See
`doc/gradient.md`.
"""
function gradient(m::Mesh, array::AbstractString; operator=:gradient,
                  method=:green_gauss, location=:cell, output::AbstractString="",
                  component::Integer=-1, overwrite::Bool=false)
    skipped = Ref{Int64}(0)
    fallback = Ref{Int64}(0)
    ptr = ccall(_sym(:mio_gradient), Ptr{Cvoid},
                (Ptr{Cvoid}, Cstring, Cstring, Cstring, Cstring, Cstring, Cint, Cint,
                 Ptr{Int64}, Ptr{Int64}),
                _handle(m), array, _op_name(operator), _method_name(method),
                String(location), output, Cint(component),
                overwrite ? Cint(1) : Cint(0), skipped, fallback)
    (mesh=Mesh(_check_ptr(ptr)), num_skipped=Int(skipped[]),
     num_fallback=Int(fallback[]))
end

# The C API spells these with hyphens; Julia symbols cannot carry one, so the
# idiomatic underscore form is translated here rather than forcing callers to
# pass strings.
_op_name(op) = String(op)
_method_name(m) = replace(String(m), '_' => '-')

"""
    hessian(mesh, array; method=:green_gauss, location=:cell, output="",
            overwrite=false) -> (; mesh, num_skipped, num_fallback)

The Hessian (second derivative) of a **scalar point-data** field —
[`gradient`](@ref)'s companion one order further, for curvature-based
adaptive refinement.

A composition of TWO `gradient` calls, not a new numerical kernel: the field
is differentiated once (point location), then that `(n, 3)` gradient is
differentiated again with the default `:gradient` operator, producing
`(n, 9)` — the flattened row-major 3x3 Hessian, `H[i,j]` at index `i*3+j`.
`method` is forwarded to BOTH internal passes. The result is named
`"<array>:hessian"` unless `output` overrides it.

A field that is at most LINEAR has an exactly zero Hessian everywhere — the
one mesh-shape-independent guarantee. For a genuinely quadratic field the
composition is exact on a structured/symmetric mesh away from its own
boundary and a good, standard, but genuinely approximate curvature estimate
on an irregular mesh (see `doc/hessian.md`). Input must have exactly one
component — a vector field's Hessian is a separate quantity per component.

Cells that cannot be evaluated yield `NaN` and are counted in `num_skipped`;
least-squares cells with a degenerate neighbourhood in either internal pass
fall back to Green-Gauss and are counted in `num_fallback` (summed over both
passes).

A curvature-driven refinement indicator needs no new function: `norm(...)`
in [`data_calc`](@ref) on the 9-component output is exactly its Frobenius
norm, ready for [`refine`](@ref)'s `where` selector.
"""
function hessian(m::Mesh, array::AbstractString; method=:green_gauss, location=:cell,
                 output::AbstractString="", overwrite::Bool=false)
    skipped = Ref{Int64}(0)
    fallback = Ref{Int64}(0)
    ptr = ccall(_sym(:mio_hessian), Ptr{Cvoid},
                (Ptr{Cvoid}, Cstring, Cstring, Cstring, Cstring, Cint, Ptr{Int64}, Ptr{Int64}),
                _handle(m), array, _method_name(method), String(location), output,
                overwrite ? Cint(1) : Cint(0), skipped, fallback)
    (mesh=Mesh(_check_ptr(ptr)), num_skipped=Int(skipped[]), num_fallback=Int(fallback[]))
end

"""
    estimate_error(mesh, array; method=:zz, marking=:none, marking_value=0.0,
                   output="", marked="", overwrite=false)
        -> (; mesh, global_error, num_skipped, num_marked)

Estimate the per-cell recovered-gradient (Zienkiewicz-Zhu) error of a
**point-data** field, and optionally mark cells for refinement.

A composition of [`gradient`](@ref) (Green-Gauss, cell location) with the
measure-weighted point<->cell averaging round trip: the indicator is
`sqrt(|measure| * sum((recovered - raw)^2))` per cell, attached as `output`
(default `"error:zz"`). `marking` is `:none` (default), `:absolute`,
`:fraction`, or `:dorfler`; when not `:none` a second Int64 0/1 array `marked`
(default `"error:marked"`) is attached too, so [`refine`](@ref)'s own `where`
selector needs no change at all — the intended use is
`refine(mesh, where="error:marked > 0.5")`. `marking_value`'s meaning depends
on `marking`: an absolute indicator threshold, a fraction in `(0, 1]` of
cells, or the Dörfler bulk fraction theta in `(0, 1]`.

Cells that cannot be evaluated read `NaN` in the indicator array and `0`
(never `NaN`) in the marking array, and are counted in `num_skipped`
(excluded from `global_error` and from `num_marked`).

Naming a *cell-data* array fails by name, for the same reason `gradient` does.
See `doc/gradient.md`.
"""
function estimate_error(m::Mesh, array::AbstractString; method=:zz, marking=:none,
                        marking_value::Real=0.0, output::AbstractString="",
                        marked::AbstractString="", overwrite::Bool=false)
    global_error = Ref{Cdouble}(0.0)
    skipped = Ref{Int64}(0)
    nmarked = Ref{Int64}(0)
    ptr = ccall(_sym(:mio_estimate_error), Ptr{Cvoid},
                (Ptr{Cvoid}, Cstring, Cstring, Cstring, Cdouble, Cstring, Cstring, Cint,
                 Ptr{Cdouble}, Ptr{Int64}, Ptr{Int64}),
                _handle(m), array, String(method), String(marking), Cdouble(marking_value),
                output, marked, overwrite ? Cint(1) : Cint(0), global_error, skipped, nmarked)
    (mesh=Mesh(_check_ptr(ptr)), global_error=Float64(global_error[]),
     num_skipped=Int(skipped[]), num_marked=Int(nmarked[]))
end

"""
    error_size(mesh; target_error, error_array="error:zz", order=1.0, min_size=0.0,
               max_size=0.0, output="size", overwrite=false)
        -> (; mesh, num_cells, num_skipped, num_clamped, min_size, max_size)

A per-point size field from a per-cell error indicator (roadmap 6.10.2), by the
adaptive procedure of Zienkiewicz and Zhu (1987): with `N` cells evaluated and
the target global error `target_error`, a cell of diameter `h` whose indicator
`error_array` (what [`estimate_error`](@ref) writes) is `eta` asks
`h * (target_error / sqrt(N) / eta)^(1 / order)`, clamped to
`[min_size, max_size]` (0: no bound); a point takes the least of its cells'
(`NaN` where none is evaluated). The sizes are written as point data `output`,
what [`mesh_cad_surface`](@ref)'s `background` reads. See `doc/error.md`.
"""
function error_size(m::Mesh; target_error::Real, error_array::AbstractString="error:zz",
                    order::Real=1.0, min_size::Real=0.0, max_size::Real=0.0,
                    output::AbstractString="size", overwrite::Bool=false)
    aname = Vector{UInt8}(codeunits(String(error_array)))
    push!(aname, 0x00)
    oname = Vector{UInt8}(codeunits(String(output)))
    push!(oname, 0x00)
    report = Ref{_CErrorSizeReport}()
    ptr = GC.@preserve aname oname begin
        opts = _CErrorSizeOpts(pointer(aname), Cdouble(target_error), Cdouble(order),
                               Cdouble(min_size), Cdouble(max_size), pointer(oname),
                               overwrite ? Int32(1) : Int32(0), Int32(0), (0, 0, 0))
        ccall(_sym(:mio_error_size), Ptr{Cvoid},
              (Ptr{Cvoid}, Ref{_CErrorSizeOpts}, Ref{_CErrorSizeReport}),
              _handle(m), Ref(opts), report)
    end
    out = Mesh(_check_ptr(ptr))
    r = report[]
    (mesh=out, num_cells=Int(r.num_cells), num_skipped=Int(r.num_skipped),
     num_clamped=Int(r.num_clamped), min_size=Float64(r.min_size),
     max_size=Float64(r.max_size))
end

const _CONFORMITY_MODES = ("conforming", "weld", "proximity")

"""
    check_conformity(mesh; tolerance=1e-9, attach_flags=true)
        -> (; mesh, dimension, num_hanging_nodes, num_t_junctions, num_partial_overlaps,
              num_unwelded_facets, num_boundary_facets, num_non_manifold_facets,
              conforming, interface_mode, warnings)

Whether `mesh` is conforming (roadmap 3.2.1): its volume cells (else its 2-D
cells) are checked for hanging nodes (a point on a volume cell's face or edge
that the cell does not use), T-junctions (the same on a 2-D cell's edge), and
partly overlapping or unwelded pairs of unmatched facets. `tolerance` is a
fraction of the examined facets' bounding-box diagonal. With `attach_flags` the
mesh carries the `conformity:hanging` point data (1 hanging node, 2 T-junction).
`interface_mode` is `"conforming"`, `"weld"` or `"proximity"`: the
[`find_interface`](@ref) mode that suits the mesh. See `doc/check_conformity.md`.
"""
function check_conformity(m::Mesh; tolerance::Real=1e-9, attach_flags::Bool=true)
    report = Ref{_CConformityReport}()
    opts = _CConformityOpts(Cdouble(tolerance), attach_flags ? Int32(1) : Int32(0), Int32(0),
                            (0, 0, 0, 0))
    ptr = ccall(_sym(:mio_check_conformity), Ptr{Cvoid},
                (Ptr{Cvoid}, Ref{_CConformityOpts}, Ref{_CConformityReport}),
                _handle(m), Ref(opts), report)
    out = Mesh(_check_ptr(ptr))
    text = _getstring((buf, n) -> ccall(_sym(:mio_check_conformity_warnings), Int64,
                                        (Ptr{UInt8}, Int64), buf, n))
    r = report[]
    (mesh=out, dimension=Int(r.dimension), num_hanging_nodes=Int(r.num_hanging_nodes),
     num_t_junctions=Int(r.num_t_junctions), num_partial_overlaps=Int(r.num_partial_overlaps),
     num_unwelded_facets=Int(r.num_unwelded_facets),
     num_boundary_facets=Int(r.num_boundary_facets),
     num_non_manifold_facets=Int(r.num_non_manifold_facets), conforming=r.conforming != 0,
     interface_mode=0 <= r.interface_mode < length(_CONFORMITY_MODES) ?
                    _CONFORMITY_MODES[r.interface_mode + 1] : string(r.interface_mode),
     warnings=isempty(text) ? String[] : String.(Base.split(text, '\n')))
end

"""
    check_tetrahedra(mesh; zero_volume_factor=1e-6)
        -> (; num_tetrahedra, num_inverted, num_degenerate, num_overlaps, min_volume,
              max_volume, valid, warnings)

Whether the tetrahedra of `mesh` are usable: its `tetra` and `tetra10` cells (a quadratic
one by its four corners) are counted as inverted (a signed volume below minus the tolerance),
degenerate (a volume within it) and, as `num_overlaps`, the pairs of tetrahedra that share a
face, used by no other, with their opposite nodes on the same side of it. The tolerance is
`zero_volume_factor` times the longest edge cubed plus the effect of six significant digits of
written coordinates. `valid` is true when there is none of the three. The mesh is not
changed. See `doc/check_tetrahedra.md`.
"""
function check_tetrahedra(m::Mesh; zero_volume_factor::Real=1e-6)
    report = Ref{_CCheckTetrahedraReport}()
    opts = _CCheckTetrahedraOpts(Cdouble(zero_volume_factor), (0, 0, 0, 0))
    _check(ccall(_sym(:mio_check_tetrahedra), Cint,
                 (Ptr{Cvoid}, Ref{_CCheckTetrahedraOpts}, Ptr{_CCheckTetrahedraReport}),
                 _handle(m), Ref(opts), report))
    text = _getstring((buf, n) -> ccall(_sym(:mio_check_tetrahedra_warnings), Int64,
                                        (Ptr{UInt8}, Int64), buf, n))
    r = report[]
    (num_tetrahedra=Int(r.num_tetrahedra), num_inverted=Int(r.num_inverted),
     num_degenerate=Int(r.num_degenerate), num_overlaps=Int(r.num_overlaps),
     min_volume=Float64(r.min_volume), max_volume=Float64(r.max_volume),
     valid=r.valid != 0,
     warnings=isempty(text) ? String[] : String.(Base.split(text, '\n')))
end

"""
    remesh(mesh, num_clusters; subdivide=-1, subsample_ratio=10.0,
           max_subdivide=4, max_iterations=100, max_repair_passes=10,
           metric=:isotropic, gradation=0.0, preserve_boundary=true,
           max_anisotropy=4.0)
        -> (; mesh, num_clusters, num_iterations, subdivide_applied,
             num_isolated_clusters, num_non_manifold_vertices)

Remesh a surface uniformly by approximated centroidal Voronoi diagram (ACVD)
clustering: replace its triangulation with a new, near-uniformly-sized,
well-shaped one at `num_clusters` vertices. Unlike every other
resolution-changing operation, the output has NEW points and NEW
connectivity with no correspondence to `mesh` — point/cell data and named
regions are dropped, field data is carried.

`metric` is `:isotropic` (default; area-weighted centroidal distance, fast,
rounds sharp features), `:quadric` (Garland-Heckbert quadric error,
preserves sharp edges/corners at extra cost per candidate move), or
`:anisotropic` (clusters shaped by a local curvature tensor — elongated
along low-curvature directions — see `max_anisotropy`). `subdivide` defaults
to automatic (`-1`): the smallest count of uniform `refine` passes reaching
`subsample_ratio` items per cluster, capped at `max_subdivide`; `0` disables
subdivision. `gradation` is the curvature-gradation exponent `gamma` in the
item weight `area * kappa^gamma`; `0.0` (default) disables gradation
entirely and reproduces plain area weighting. `preserve_boundary` (default
`true`) detects the input's open boundary (if any), seeds it before the
interior, and emits a `line` dual cell along boundary edges whose endpoints
land in different clusters — a no-op on a closed mesh. `max_anisotropy`
(default `4.0`, a measured value — see `kRemeshDefaultMaxAnisotropy`'s doc
comment in `remesh.hpp`) is, under `metric=:anisotropic`, the maximum ratio
between the two in-plane target edge lengths a curvature tensor may
request; `1.0` recovers the isotropic shape exactly. Must be at least `1.0`,
and an error to set away from the default under any other metric.
`num_isolated_clusters` and `num_non_manifold_vertices` are two distinct
causes of non-manifold output (disconnected clusters vs. "bowtie" vertices)
that repair could not fully fix; check both rather than assuming.

Goes through `mio_remesh_ex`/`mio_remesh_opts` rather than the flat
`mio_remesh` — the `refine`/`mio_refine_ex` precedent, needed because
`mio_remesh` is a flat C function with no room to grow (it already changed
once).

Attribution: the isotropic clustering engine is derived from
[pyacvd](https://github.com/pyvista/pyacvd) (MIT, (c) 2017-2024 The PyVista
Developers), itself an independent implementation of the published research
of S. Valette and J.-M. Chassery, not of the CeCILL-B licensed ACVD project.
See `doc/remesh.md`.
"""
function remesh(m::Mesh, num_clusters::Integer; subdivide::Integer=-1,
                subsample_ratio::Real=10.0, max_subdivide::Integer=4,
                max_iterations::Integer=100, max_repair_passes::Integer=10,
                metric::Symbol=:isotropic, gradation::Real=0.0,
                preserve_boundary::Bool=true, max_anisotropy::Real=4.0)
    metric_c = Vector{UInt8}(codeunits(String(metric)))
    push!(metric_c, 0x00)
    report = Ref{_CRemeshReport}()
    ptr = GC.@preserve metric_c begin
        opts = _CRemeshOpts(Int64(num_clusters), Int32(subdivide), Int32(max_subdivide),
                            Cdouble(subsample_ratio), Int32(max_iterations),
                            Int32(max_repair_passes), Cstring(pointer(metric_c)),
                            Cdouble(gradation), preserve_boundary ? Int32(1) : Int32(0),
                            Int32(0), Cdouble(max_anisotropy),
                            (Cdouble(0), Cdouble(0), Cdouble(0)),
                            (Int64(0), Int64(0), Int64(0), Int64(0)))
        ccall(_sym(:mio_remesh_ex), Ptr{Cvoid},
              (Ptr{Cvoid}, Ref{_CRemeshOpts}, Ptr{_CRemeshReport}),
              _handle(m), Ref(opts), report)
    end
    r = _check_ptr(ptr)
    rep = report[]
    (mesh=Mesh(r), num_clusters=Int(rep.num_clusters),
     num_iterations=Int(rep.num_iterations), subdivide_applied=Int(rep.subdivide_applied),
     num_isolated_clusters=Int(rep.num_isolated_clusters),
     num_non_manifold_vertices=Int(rep.num_non_manifold_vertices))
end

"""
    remesh_volume(mesh; resolution=nothing, cell_size=nothing, bounds=nothing,
                  padding=0.0, padding_relative=0.1, max_cells=20_000_000,
                  max_tets=20_000_000, warp_fraction=0.35, sign=:pseudonormal,
                  watertight_check=:warn)
        -> (; mesh, num_tets, num_vertices_warped, num_tets_rejected,
             num_non_manifold_edges, input_quality)

Retetrahedralize `mesh` (a volume mesh, or a closed surface) at a
caller-chosen resolution by isosurface stuffing — the volumetric sibling of
[`remesh`](@ref). Unlike `remesh`, `mesh` may be a VOLUME mesh directly (its
boundary is extracted internally) as well as a closed surface. Same
no-correspondence-with-the-input output contract as `remesh`: new points, new
connectivity, `point_data`/`cell_data`/named regions dropped, `field_data`
carried.

Exactly one of `resolution`/`cell_size` must be given, sizing a
body-centered cubic (BCC) lattice whose uncut tets have dihedral angles from
a fixed, mesh-size-independent set. `warp_fraction` (default `0.35`) moves
lattice vertices near the surface onto it, trading a small, measured chance
of non-manifold boundary edges (reported as `num_non_manifold_edges`) for
substantially better boundary tet quality; `0` disables warping and gives an
exactly watertight but lower-quality boundary. See `doc/remesh_volume.md`
for the measured tradeoff.

Attribution: implemented from the PUBLISHED DESCRIPTION of Labelle &
Shewchuk, "Isosurface Stuffing" (SIGGRAPH 2007) only — neither the paper's
own reference implementation (Stellar) nor TetGen (AGPL-3.0) is read or
vendored here. See `doc/remesh_volume.md`.
"""
function remesh_volume(m::Mesh; resolution=nothing, cell_size=nothing, bounds=nothing,
                       padding::Real=0.0, padding_relative::Real=0.1,
                       max_cells::Integer=20_000_000, max_tets::Integer=20_000_000,
                       warp_fraction::Real=0.35, sign=:pseudonormal, weight=:angle,
                       watertight_check=:warn, surface_region::AbstractString="",
                       grid_cell_size::Real=0.0, max_winding_work::Real=2.0e9)
    res = resolution === nothing ? Int64[] : Int64[Int64(v) for v in resolution]
    (resolution === nothing || length(res) == 3) ||
        throw(ArgumentError("resolution must have three entries"))
    bnd = bounds === nothing ? Cdouble[] : Cdouble[Cdouble(v) for v in bounds]
    (bounds === nothing || length(bnd) == 6) ||
        throw(ArgumentError("bounds must have six entries"))

    region_buf = Vector{UInt8}(codeunits(String(surface_region)))
    push!(region_buf, 0x00)
    report = Ref{_CRemeshVolumeReport}()
    ptr = GC.@preserve res bnd region_buf begin
        inner = _sdf_opts(; sign, weight, location=:corner, band=0.0,
                          record_closest_cell=false, record_inside=false, watertight_check,
                          surface_region, grid_cell_size, max_winding_work, region_buf)
        opts = _CRemeshVolumeOpts(isempty(res) ? Ptr{Int64}(C_NULL) : pointer(res),
                                  isempty(bnd) ? Ptr{Cdouble}(C_NULL) : pointer(bnd),
                                  Cdouble(cell_size === nothing ? 0.0 : cell_size),
                                  Cdouble(padding), Cdouble(padding_relative),
                                  Int64(max_cells), Int64(max_tets), Cdouble(warp_fraction),
                                  (Int64(0), Int64(0), Int64(0), Int64(0), Int64(0), Int64(0)),
                                  inner)
        ccall(_sym(:mio_remesh_volume_ex), Ptr{Cvoid},
              (Ptr{Cvoid}, Ref{_CRemeshVolumeOpts}, Ptr{_CRemeshVolumeReport}),
              _handle(m), Ref(opts), report)
    end
    r = _check_ptr(ptr)
    rep = report[]
    q = rep.input_quality
    (mesh=Mesh(r), num_tets=Int(rep.num_tets),
     num_vertices_warped=Int(rep.num_vertices_warped),
     num_tets_rejected=Int(rep.num_tets_rejected),
     num_non_manifold_edges=Int(rep.num_non_manifold_edges),
     input_quality=(boundary_edges=Int(q.boundary_edges),
                    non_manifold_edges=Int(q.non_manifold_edges),
                    inconsistent_pairs=Int(q.inconsistent_pairs),
                    degenerate_triangles=Int(q.degenerate_triangles),
                    watertight=q.watertight != 0))
end

# --- combining / comparing ---------------------------------------------------

"""
    merge(meshes; weld=false, atol=1e-12, source_tag=true,
          data_policy=:intersection, drop_duplicate_cells=false) -> Mesh

Combine two or more meshes into one. `data_policy` is `:intersection` (keep
only keys present in every input) or `:fill` (union, NaN for missing rows).
`source_tag` adds an Int64 `"source_mesh_id"` cell-data array.

Not exported (it would shadow `Base.merge`): call it qualified.
"""
function merge(meshes::AbstractVector{Mesh}; weld::Bool=false, atol::Real=1e-12,
               source_tag::Bool=true, data_policy::Symbol=:intersection,
               drop_duplicate_cells::Bool=false)
    isempty(meshes) && throw(ArgumentError("merge needs at least one mesh"))
    policy = data_policy === :intersection ? Cint(0) : data_policy === :fill ? Cint(1) :
             throw(ArgumentError("data_policy must be :intersection or :fill"))
    handles = Ptr{Cvoid}[_handle(x) for x in meshes]
    ptr = GC.@preserve meshes handles ccall(
        _sym(:mio_merge), Ptr{Cvoid},
        (Ptr{Ptr{Cvoid}}, Int64, Cint, Cdouble, Cint, Cint, Cint),
        pointer(handles), Int64(length(handles)), weld ? Cint(1) : Cint(0),
        Float64(atol), source_tag ? Cint(1) : Cint(0), policy,
        drop_duplicate_cells ? Cint(1) : Cint(0))
    Mesh(_check_ptr(ptr))
end

"""
    undo_green(coarse, fine) -> (; mesh, num_groups_undone, num_cells_removed)

Green-element undo: restore `fine`'s transitional (closure-only) cells back to
their original parent, read verbatim from `coarse` — a lookup, not a
reconstruction, since `refine` never renumbers or prunes points. Undoes
`refine`'s known quality-degradation issue with repeated selective passes over
the same region ("restore the parent and re-split from scratch").

`fine` must carry `refine:cell_id`/`refine:parent_id`/`refine:level` (i.e. it
must come from `refine(coarse, ...; record_hierarchy=true,
record_levels=true)`); `coarse` must be the exact mesh that call was run on.
The six reserved `refine:*` arrays are dropped from the output. Only a
single-pass (`levels=1`) hierarchy is supported. See `doc/undo_green.md`.
"""
function undo_green(coarse::Mesh, fine::Mesh)
    ngroups = Ref{Int64}(0); nremoved = Ref{Int64}(0)
    ptr = ccall(_sym(:mio_undo_green), Ptr{Cvoid},
                (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Int64}, Ptr{Int64}),
                _handle(coarse), _handle(fine), ngroups, nremoved)
    (mesh=Mesh(_check_ptr(ptr)), num_groups_undone=Int(ngroups[]),
     num_cells_removed=Int(nremoved[]))
end

"""
    boundary_layer(mesh; first_height, num_layers=3, growth=1.2, reverse=false,
                   collision_fraction=0.45, smooth_iterations=5, max_shrink=6,
                   collapse="points", normal_smooth_iterations=0, front_check=true)
        -> (; mesh, num_cells, num_wedges, num_hexahedra, num_pyramids, num_tetrahedra,
              num_shrunk, num_collapsed, num_not_visible, num_fewer_layers,
              num_front_crossings, min_scaled_jacobian, min_thickness, max_thickness,
              cap_closed)

Grow prismatic boundary layers from a triangle and quad surface (roadmap 6.9):
wedges from triangles and hexahedra from quads along each point's most normal
normal (Aubry and Löhner 2009), the first `first_height` thick and each next
`growth` times the one before; stopped at `collision_fraction` of the distance
to the surface ahead; the points of a cell that would fold shrunk and then
collapsed into pyramids and tetrahedra (`collapse="points"`) or their top layer
dropped (`collapse="layers"`, counted by `num_fewer_layers`).
`normal_smooth_iterations` Laplacian passes smooth the directions, and
`front_check` tests the last layer's faces against each other and the surface
(`num_front_crossings`). The wall and the cap come out as faces, with
`boundary_layer:layer` cell data. See `doc/boundary_layer.md`.
"""
function boundary_layer(mesh::Mesh; first_height::Real, num_layers::Integer=3,
                        growth::Real=1.2, reverse::Bool=false, collision_fraction::Real=0.45,
                        smooth_iterations::Integer=5, max_shrink::Integer=6,
                        collapse::AbstractString="points",
                        normal_smooth_iterations::Integer=0, front_check::Bool=true)
    collapse_c = collapse == "points" ? Int32(0) : collapse == "layers" ? Int32(1) :
                 throw(ArgumentError("meshio++: boundary_layer: collapse must be \"points\" or \"layers\""))
    opts = _CBoundaryLayerOpts(Cdouble(first_height), Cdouble(growth),
                               Cdouble(collision_fraction), Int32(num_layers), Int32(reverse),
                               Int32(smooth_iterations), Int32(max_shrink), collapse_c,
                               Int32(normal_smooth_iterations), Int32(front_check), Int32(0),
                               (0,))
    report = Ref{_CBoundaryLayerReport}()
    ptr = ccall(_sym(:mio_boundary_layer), Ptr{Cvoid},
                (Ptr{Cvoid}, Ref{_CBoundaryLayerOpts}, Ref{_CBoundaryLayerReport}),
                _handle(mesh), Ref(opts), report)
    out = Mesh(_check_ptr(ptr))
    r = report[]
    (mesh=out, num_cells=Int(r.num_cells), num_wedges=Int(r.num_wedges),
     num_hexahedra=Int(r.num_hexahedra), num_pyramids=Int(r.num_pyramids),
     num_tetrahedra=Int(r.num_tetrahedra), num_shrunk=Int(r.num_shrunk),
     num_collapsed=Int(r.num_collapsed), num_not_visible=Int(r.num_not_visible),
     num_fewer_layers=Int(r.num_fewer_layers), num_front_crossings=Int(r.num_front_crossings),
     min_scaled_jacobian=Float64(r.min_scaled_jacobian),
     min_thickness=Float64(r.min_thickness), max_thickness=Float64(r.max_thickness),
     cap_closed=r.cap_closed != 0)
end

_vec3(v, name) = length(v) == 3 ? (Cdouble(v[1]), Cdouble(v[2]), Cdouble(v[3])) :
                 throw(ArgumentError("meshio++: $name must have 3 entries"))

function _sweep_result(out::Mesh, r::_CSweepReport, text::AbstractString)
    (mesh=out, num_cells=Int(r.num_cells), num_hexahedra=Int(r.num_hexahedra),
     num_wedges=Int(r.num_wedges), num_pyramids=Int(r.num_pyramids),
     num_tetrahedra=Int(r.num_tetrahedra), num_quads=Int(r.num_quads),
     num_triangles=Int(r.num_triangles), num_flipped=Int(r.num_flipped),
     num_on_axis=Int(r.num_on_axis), min_scaled_jacobian=Float64(r.min_scaled_jacobian),
     full=r.full != 0, warnings=isempty(text) ? String[] : String.(Base.split(text, '\n')))
end

"""
    extrude(mesh; direction=(0.0, 0.0, 1.0), num_layers=1, ratio=1.0, offsets=nothing,
            record_layers=false, max_cells=20_000_000)
        -> (; mesh, num_cells, num_hexahedra, num_wedges, num_pyramids, num_tetrahedra,
              num_quads, num_triangles, num_flipped, num_on_axis, min_scaled_jacobian,
              full, warnings)

Sweep a surface's triangles and quads into wedges and hexahedra (or, with no
2-D cells, its lines into quads) along `direction`, the whole displacement, in
`num_layers` layers spaced uniformly, by the geometric `ratio`, or at the
explicit `offsets` (fractions, increasing, ending at 1) (roadmap 3.1.2). The
result holds the swept cells only, with `Side` regions `bottom`, `top` and
`side`. See `doc/extrude.md`.
"""
function extrude(m::Mesh; direction=(0.0, 0.0, 1.0), num_layers::Integer=1, ratio::Real=1.0,
                 offsets=nothing, record_layers::Bool=false, max_cells::Integer=20_000_000)
    offs = offsets === nothing ? Cdouble[] : Cdouble.(collect(offsets))
    report = Ref{_CSweepReport}()
    ptr = GC.@preserve offs begin
        opts = _CExtrudeOpts(_vec3(direction, "extrude: direction"), Cdouble(ratio),
                             isempty(offs) ? Ptr{Cdouble}(C_NULL) : pointer(offs),
                             Int64(length(offs)), Int64(max_cells), Int32(num_layers),
                             Int32(record_layers), (0, 0, 0, 0))
        ccall(_sym(:mio_extrude), Ptr{Cvoid},
              (Ptr{Cvoid}, Ref{_CExtrudeOpts}, Ref{_CSweepReport}),
              _handle(m), Ref(opts), report)
    end
    out = Mesh(_check_ptr(ptr))
    text = _getstring((buf, n) -> ccall(_sym(:mio_extrude_warnings), Int64,
                                        (Ptr{UInt8}, Int64), buf, n))
    _sweep_result(out, report[], text)
end

"""
    revolve(mesh; axis=(0.0, 1.0, 0.0), axis_point=(0.0, 0.0, 0.0), angle=360.0,
            num_layers=16, tolerance=1e-9, record_layers=false, max_cells=20_000_000)
        -> the named tuple of [`extrude`](@ref)

Sweep a surface (or a curve) around the axis through `axis_point` along `axis`
by `angle` degrees (right-hand rule; 360 closes the turn) in `num_layers` steps
(roadmap 3.1.3). Points within `tolerance` (a fraction of the bounding-box
diagonal) of the axis stay put, so the cells around them collapse into
pyramids, tetrahedra and wedges. `Side` regions `start`, `end` (a partial turn)
and `side`. See `doc/extrude.md`.
"""
function revolve(m::Mesh; axis=(0.0, 1.0, 0.0), axis_point=(0.0, 0.0, 0.0), angle::Real=360.0,
                 num_layers::Integer=16, tolerance::Real=1e-9, record_layers::Bool=false,
                 max_cells::Integer=20_000_000)
    opts = _CRevolveOpts(_vec3(axis_point, "revolve: axis_point"), _vec3(axis, "revolve: axis"),
                         Cdouble(angle), Cdouble(tolerance), Int64(max_cells),
                         Int32(num_layers), Int32(record_layers), (0, 0, 0, 0))
    report = Ref{_CSweepReport}()
    ptr = ccall(_sym(:mio_revolve), Ptr{Cvoid},
                (Ptr{Cvoid}, Ref{_CRevolveOpts}, Ref{_CSweepReport}),
                _handle(m), Ref(opts), report)
    out = Mesh(_check_ptr(ptr))
    text = _getstring((buf, n) -> ccall(_sym(:mio_revolve_warnings), Int64,
                                        (Ptr{UInt8}, Int64), buf, n))
    _sweep_result(out, report[], text)
end

"""
    mirror(mesh; origin=(0.0, 0.0, 0.0), normal=(1.0, 0.0, 0.0), keep_original=false,
           weld=true, tolerance=1e-9, reflect_vector_data=false)
        -> (; mesh, num_flipped_cells, num_points_on_plane, num_dropped_cells,
              num_dropped_sides, warnings)

Reflect `mesh` across the plane through `origin` with `normal`, every cell
right side out (roadmap 3.2.5): nodes permuted, `Side` region entries moved to
their facets' new numbers. With `keep_original` the input followed by its
image, the points on the plane shared when `weld`; `reflect_vector_data`
reflects float vector and tensor data. See `doc/mirror.md`.
"""
function mirror(m::Mesh; origin=(0.0, 0.0, 0.0), normal=(1.0, 0.0, 0.0),
                keep_original::Bool=false, weld::Bool=true, tolerance::Real=1e-9,
                reflect_vector_data::Bool=false)
    opts = _CMirrorOpts(_vec3(origin, "mirror: origin"), _vec3(normal, "mirror: normal"),
                        Cdouble(tolerance), Int32(keep_original), Int32(weld),
                        Int32(reflect_vector_data), Int32(0), (0, 0, 0, 0))
    report = Ref{_CMirrorReport}()
    ptr = ccall(_sym(:mio_mirror), Ptr{Cvoid},
                (Ptr{Cvoid}, Ref{_CMirrorOpts}, Ref{_CMirrorReport}),
                _handle(m), Ref(opts), report)
    out = Mesh(_check_ptr(ptr))
    text = _getstring((buf, n) -> ccall(_sym(:mio_mirror_warnings), Int64,
                                        (Ptr{UInt8}, Int64), buf, n))
    r = report[]
    (mesh=out, num_flipped_cells=Int(r.num_flipped_cells),
     num_points_on_plane=Int(r.num_points_on_plane), num_dropped_cells=Int(r.num_dropped_cells),
     num_dropped_sides=Int(r.num_dropped_sides),
     warnings=isempty(text) ? String[] : String.(Base.split(text, '\n')))
end

"""
    find_intersections(mesh, other=nothing; region="", region_other="", segments=true,
                       max_pairs=0)
        -> (; mesh, other, segments, pairs, num_pairs, num_triangle_pairs, num_coplanar,
              num_segments, num_cells_a, num_cells_b, truncated, clean, warnings)

Where the surface `mesh` intersects itself, or where it and the surface
`other` interpenetrate or touch, by exact triangle-triangle tests (roadmap
3.1.2). `mesh` (and `other`) come back with the `intersections:count` cell
data, `segments` is the line mesh of the intersection segments (`nothing` when
not asked for) and `pairs` a `2 x n` matrix of **1-based** global cell ids,
at most `max_pairs` of them (0: all). See `doc/find_intersections.md`.
"""
function find_intersections(m::Mesh, other::Union{Mesh,Nothing}=nothing;
                            region::AbstractString="", region_other::AbstractString="",
                            segments::Bool=true, max_pairs::Integer=0)
    report = Ref{_CFindIntersectionsReport}()
    ra = Vector{UInt8}(codeunits(String(region) * "\0"))
    rb = Vector{UInt8}(codeunits(String(region_other) * "\0"))
    res = GC.@preserve ra rb begin
        opts = _CFindIntersectionsOpts(Cstring(pointer(ra)), Cstring(pointer(rb)),
                                       Int64(max_pairs), Int32(1), Int32(segments),
                                       (0, 0, 0, 0))
        ccall(_sym(:mio_find_intersections), Ptr{Cvoid},
              (Ptr{Cvoid}, Ptr{Cvoid}, Ref{_CFindIntersectionsOpts},
               Ref{_CFindIntersectionsReport}),
              _handle(m), other === nothing ? C_NULL : _handle(other), Ref(opts), report)
    end
    _check_ptr(res)
    local out, out_b, segs, pairs
    try
        out = Mesh(_check_ptr(ccall(_sym(:mio_find_intersections_result_take_mesh), Ptr{Cvoid},
                                    (Ptr{Cvoid},), res)))
        pb = ccall(_sym(:mio_find_intersections_result_take_mesh_b), Ptr{Cvoid},
                   (Ptr{Cvoid},), res)
        out_b = pb == C_NULL ? nothing : Mesh(pb)
        ps = ccall(_sym(:mio_find_intersections_result_take_segments), Ptr{Cvoid},
                   (Ptr{Cvoid},), res)
        segs = ps == C_NULL ? nothing : Mesh(ps)
        count = Ref{Int64}(0)
        p = ccall(_sym(:mio_find_intersections_result_pairs), Ptr{Int64},
                  (Ptr{Cvoid}, Ref{Int64}), res, count)
        pairs = Matrix{Int}(undef, 2, count[])
        for k in 1:(2 * count[])
            pairs[k] = Int(unsafe_load(p, k)) + 1
        end
    finally
        ccall(_sym(:mio_find_intersections_result_free), Cvoid, (Ptr{Cvoid},), res)
    end
    text = _getstring((buf, n) -> ccall(_sym(:mio_find_intersections_warnings), Int64,
                                        (Ptr{UInt8}, Int64), buf, n))
    r = report[]
    (mesh=out, other=out_b, segments=segs, pairs=pairs, num_pairs=Int(r.num_pairs),
     num_triangle_pairs=Int(r.num_triangle_pairs), num_coplanar=Int(r.num_coplanar),
     num_segments=Int(r.num_segments), num_cells_a=Int(r.num_cells_a),
     num_cells_b=Int(r.num_cells_b), truncated=r.truncated != 0, clean=r.clean != 0,
     warnings=isempty(text) ? String[] : String.(Base.split(text, '\n')))
end

const _THICKNESS_DIRECTIONS = Dict("auto" => Int32(0), "inward" => Int32(1),
                                   "outward" => Int32(2))
const _THICKNESS_METHODS = Dict("ray" => Int32(0), "sphere" => Int32(1), "cone" => Int32(2))

"""
    thickness(mesh, other=nothing; direction="auto", max_distance=0.0, threshold=0.0,
              region="", region_other="", method="ray", cone_angle=15.0, cone_rays=16,
              two_sided=false, max_normal_deviation=90.0)
        -> (; mesh, direction, method, num_cells, num_hit, num_missed, num_below, min, max,
              mean, p2, p95, min_cell, min_point, max_distance, warnings)

Wall thickness or clearance per facet of the surface `mesh`, by a ray from each
cell's centroid along its normal (roadmap 3.1.4): `direction` `"inward"` stops
at the far side of the wall, `"outward"` at the next wall facing the ray, and
`"auto"` is inward for one surface, outward against `other`. `mesh` comes back
with the `thickness:distance` (NaN for a miss) and `thickness:hit_cell` cell
data; `min_cell` is **1-based** (0: no hit). `method` `"sphere"` or `"cone"`,
`two_sided` and `max_normal_deviation` are roadmap 7.5's methods and scan-safe
hits. See `doc/thickness.md`.
"""
function thickness(m::Mesh, other::Union{Mesh,Nothing}=nothing;
                   direction::AbstractString="auto", max_distance::Real=0.0,
                   threshold::Real=0.0, region::AbstractString="",
                   region_other::AbstractString="", method::AbstractString="ray",
                   cone_angle::Real=15.0, cone_rays::Integer=16, two_sided::Bool=false,
                   max_normal_deviation::Real=90.0)
    haskey(_THICKNESS_DIRECTIONS, direction) ||
        throw(ArgumentError("meshio++: thickness: direction must be auto, inward or outward"))
    haskey(_THICKNESS_METHODS, method) ||
        throw(ArgumentError("meshio++: thickness: method must be ray, sphere or cone"))
    report = Ref{_CThicknessReport}()
    ra = Vector{UInt8}(codeunits(String(region) * "\0"))
    rb = Vector{UInt8}(codeunits(String(region_other) * "\0"))
    ptr = GC.@preserve ra rb begin
        opts = _CThicknessOpts(Cstring(pointer(ra)), Cstring(pointer(rb)), Cdouble(max_distance),
                               Cdouble(threshold), _THICKNESS_DIRECTIONS[direction], Int32(1),
                               _THICKNESS_METHODS[method], two_sided ? Int32(1) : Int32(0),
                               Cdouble(cone_angle), Int64(cone_rays),
                               Cdouble(max_normal_deviation))
        ccall(_sym(:mio_thickness), Ptr{Cvoid},
              (Ptr{Cvoid}, Ptr{Cvoid}, Ref{_CThicknessOpts}, Ref{_CThicknessReport}),
              _handle(m), other === nothing ? C_NULL : _handle(other), Ref(opts), report)
    end
    out = Mesh(_check_ptr(ptr))
    text = _getstring((buf, n) -> ccall(_sym(:mio_thickness_warnings), Int64,
                                        (Ptr{UInt8}, Int64), buf, n))
    r = report[]
    (mesh=out, direction=r.direction == 2 ? "outward" : "inward",
     method=("ray", "sphere", "cone")[r.method + 1], num_cells=Int(r.num_cells),
     num_hit=Int(r.num_hit), num_missed=Int(r.num_missed), num_below=Int(r.num_below),
     min=r.min, max=r.max, mean=r.mean, p2=r.p2, p95=r.p95, min_cell=r.min_cell < 0 ? 0 : Int(r.min_cell) + 1,
     min_point=collect(r.min_point), max_distance=r.max_distance,
     warnings=isempty(text) ? String[] : String.(Base.split(text, '\n')))
end

const _DEVIATION_TYPES = Dict("plane" => Int32(1), "sphere" => Int32(2), "cylinder" => Int32(3),
                               "cone" => Int32(4), "torus" => Int32(5))

function _deviation_result(ptr, report)
    out = Mesh(_check_ptr(ptr))
    text = _getstring((buf, n) -> ccall(_sym(:mio_deviation_warnings), Int64,
                                        (Ptr{UInt8}, Int64), buf, n))
    r = report[]
    nb = Int(r.num_bins)
    counts = zeros(Int64, nb)
    edges = zeros(Cdouble, nb + 1)
    ccall(_sym(:mio_deviation_histogram), Int64, (Ptr{Int64}, Ptr{Cdouble}, Int64),
          counts, edges, nb)
    (mesh=out, num_points=Int(r.num_points), num_measured=Int(r.num_measured),
     min=r.min, max=r.max, mean=r.mean, rms=r.rms, sigma=r.sigma,
     percentiles=(p1=r.percentiles[1], p5=r.percentiles[2], p50=r.percentiles[3],
                  p95=r.percentiles[4], p99=r.percentiles[5]),
     histogram_edges=edges, histogram_counts=counts,
     num_outside_range=Int(r.num_outside_range), share_within=r.share_within,
     suggested_range=r.suggested_range, num_without_normal=Int(r.num_without_normal),
     warnings=isempty(text) ? String[] : String.(Base.split(text, '\n')),
     registered=r.registered != 0, registration_rms=r.registration_rms,
     matrix=_deviation_matrix())
end

function _deviation_matrix()
    m = zeros(Cdouble, 16)
    ccall(_sym(:mio_deviation_matrix), Int64, (Ptr{Cdouble},), m)
    [m[4 * (i - 1) + j] for i in 1:4, j in 1:4]
end

"""
    deviation(scan, reference; tolerance=0.0, bins=20, range=nothing,
              max_normal_angle=90.0, reference_region="", register=false)
    deviation(scan, type, params; tolerance=0.0, bins=20, range=nothing)
        -> (; mesh, num_points, num_measured, min, max, mean, rms, sigma, percentiles,
              histogram_edges, histogram_counts, num_outside_range, share_within,
              suggested_range, num_without_normal, warnings, registered,
              registration_rms, matrix)

The signed deviation of every point of `scan` from the surface `reference`
(positive outside it), or from a fitted primitive `type` (`"plane"`,
`"sphere"`, `"cylinder"`, `"cone"`, `"torus"`) with its ten `fit:patch_params`
entries (roadmap 7.4). `mesh` comes back with the `deviation:distance` point
data. `max_normal_angle` is the far-side rule (180 turns it off). `register`
registers the scan first (`register_rigid` with `initial="auto"`); `mesh` is
then the registered scan and `matrix` the motion (roadmap 7.4.1). See
`doc/deviation.md`.
"""
function deviation(scan::Mesh, reference::Mesh; tolerance::Real=0.0, bins::Integer=20,
                   range=nothing, max_normal_angle::Real=90.0,
                   reference_region::AbstractString="", register::Bool=false)
    lo, hi = range === nothing ? (0.0, 0.0) : (Float64(range[1]), Float64(range[2]))
    report = Ref{_CDeviationReport}()
    rr = Vector{UInt8}(codeunits(String(reference_region) * "\0"))
    ptr = GC.@preserve rr begin
        opts = _CDeviationOpts(Cdouble(tolerance), Int64(bins), lo, hi,
                               Cdouble(max_normal_angle), Cstring(pointer(rr)), Int32(1),
                               register ? Int32(1) : Int32(0), (0, 0, 0, 0))
        ccall(_sym(:mio_deviation), Ptr{Cvoid},
              (Ptr{Cvoid}, Ptr{Cvoid}, Ref{_CDeviationOpts}, Ref{_CDeviationReport}),
              _handle(scan), _handle(reference), Ref(opts), report)
    end
    _deviation_result(ptr, report)
end

function deviation(scan::Mesh, type::AbstractString, params::AbstractVector{<:Real};
                   tolerance::Real=0.0, bins::Integer=20, range=nothing)
    haskey(_DEVIATION_TYPES, type) ||
        throw(ArgumentError("meshio++: deviation: unknown primitive type $(repr(type))"))
    length(params) == 10 ||
        throw(ArgumentError("meshio++: deviation: a primitive takes ten parameters"))
    lo, hi = range === nothing ? (0.0, 0.0) : (Float64(range[1]), Float64(range[2]))
    report = Ref{_CDeviationReport}()
    p = Cdouble.(collect(params))
    opts = _CDeviationOpts(Cdouble(tolerance), Int64(bins), lo, hi, 90.0, Cstring(C_NULL),
                           Int32(1), Int32(0), (0, 0, 0, 0))
    ptr = ccall(_sym(:mio_deviation_primitive), Ptr{Cvoid},
                (Ptr{Cvoid}, Int32, Ptr{Cdouble}, Ref{_CDeviationOpts}, Ref{_CDeviationReport}),
                _handle(scan), _DEVIATION_TYPES[type], p, Ref(opts), report)
    _deviation_result(ptr, report)
end

"""
    deviation_regions(mesh, tolerance; field="deviation:distance", prefix="deviation:",
                      min_area=0.0)
        -> (; mesh, regions, num_excess, num_missing, num_dropped, warnings)

The connected regions of a measured surface beyond `+-tolerance` (roadmap
7.4.3): excess and missing material, written as named cell regions with the
`deviation:region` cell data (0-based region indices, -1 within tolerance).
Each of `regions` (largest `|volume|` first) is `(; name, sign, num_cells,
area, extreme, mean, volume, centroid, bounds)`. See `doc/deviation.md`.
"""
function deviation_regions(m::Mesh, tolerance::Real; field::AbstractString="deviation:distance",
                           prefix::AbstractString="deviation:", min_area::Real=0.0)
    report = Ref{_CDeviationRegionsReport}()
    fb = Vector{UInt8}(codeunits(String(field) * "\0"))
    pb = Vector{UInt8}(codeunits(String(prefix) * "\0"))
    ptr = GC.@preserve fb pb begin
        opts = _CDeviationRegionsOpts(Cdouble(tolerance), Cdouble(min_area),
                                      Cstring(pointer(fb)), Cstring(pointer(pb)), (0, 0, 0, 0))
        ccall(_sym(:mio_deviation_regions), Ptr{Cvoid},
              (Ptr{Cvoid}, Ref{_CDeviationRegionsOpts}, Ref{_CDeviationRegionsReport}),
              _handle(m), Ref(opts), report)
    end
    out = Mesh(_check_ptr(ptr))
    r = report[]
    regions = map(0:(Int(r.num_regions) - 1)) do k
        info = Ref{_CDeviationRegionInfo}()
        ccall(_sym(:mio_deviation_region_get), Cint, (Int64, Ref{_CDeviationRegionInfo}), k, info)
        name = _getstring((buf, n) -> ccall(_sym(:mio_deviation_region_name), Int64,
                                            (Int64, Ptr{UInt8}, Int64), k, buf, n))
        g = info[]
        (name=name, sign=Int(g.sign), num_cells=Int(g.num_cells), area=g.area,
         extreme=g.extreme, mean=g.mean, volume=g.volume, centroid=collect(g.centroid),
         bounds=collect(g.bounds))
    end
    text = _getstring((buf, n) -> ccall(_sym(:mio_deviation_regions_warnings), Int64,
                                        (Ptr{UInt8}, Int64), buf, n))
    (mesh=out, regions=regions, num_excess=Int(r.num_excess), num_missing=Int(r.num_missing),
     num_dropped=Int(r.num_dropped),
     warnings=isempty(text) ? String[] : String.(Base.split(text, '\n')))
end

"""
    exaggerate(mesh; field="deviation:distance", scale=10.0, clamp=0.0, smooth_fraction=0.01)
        -> (; mesh, clamp, radius, num_moved, warnings)

A copy of the surface `mesh` with each point displaced along its normal by
`scale` times its clamped, smoothed deviation (roadmap 7.4.4), for a figure.
See `doc/deviation.md`.
"""
function exaggerate(m::Mesh; field::AbstractString="deviation:distance", scale::Real=10.0,
                    clamp::Real=0.0, smooth_fraction::Real=0.01)
    report = Ref{_CExaggerateReport}()
    fb = Vector{UInt8}(codeunits(String(field) * "\0"))
    ptr = GC.@preserve fb begin
        opts = _CExaggerateOpts(Cstring(pointer(fb)), Cdouble(scale), Cdouble(clamp),
                                Cdouble(smooth_fraction), (0, 0, 0, 0))
        ccall(_sym(:mio_exaggerate), Ptr{Cvoid},
              (Ptr{Cvoid}, Ref{_CExaggerateOpts}, Ref{_CExaggerateReport}),
              _handle(m), Ref(opts), report)
    end
    out = Mesh(_check_ptr(ptr))
    text = _getstring((buf, n) -> ccall(_sym(:mio_exaggerate_warnings), Int64,
                                        (Ptr{UInt8}, Int64), buf, n))
    r = report[]
    (mesh=out, clamp=r.clamp, radius=r.radius, num_moved=Int(r.num_moved),
     warnings=isempty(text) ? String[] : String.(Base.split(text, '\n')))
end

"""
    colorize(mesh, field; cmap="viridis", vmin=nothing, vmax=nothing, symmetric=false,
             steps=0, nan_color="#808080", component=nothing) -> (; mesh, vmin, vmax, location)

Bake `field` into the UInt8 data `red`, `green` and `blue` through a colormap
(`steps` of 2 or more bands it at fixed steps; roadmap 7.4.4), on the field's
location (`"point"` or `"cell"`). See `doc/colorize.md`.
"""
function colorize(m::Mesh, field::AbstractString; cmap::AbstractString="viridis",
                  vmin=nothing, vmax=nothing, symmetric::Bool=false, steps::Integer=0,
                  nan_color::AbstractString="#808080", component=nothing)
    (vmin === nothing) == (vmax === nothing) ||
        throw(ArgumentError("meshio++: colorize: give both vmin and vmax, or neither"))
    lo, hi = vmin === nothing ? (0.0, 0.0) : (Float64(vmin), Float64(vmax))
    report = Ref{_CColorizeReport}()
    fb = Vector{UInt8}(codeunits(String(field) * "\0"))
    cb = Vector{UInt8}(codeunits(String(cmap) * "\0"))
    nb = Vector{UInt8}(codeunits(String(nan_color) * "\0"))
    ptr = GC.@preserve fb cb nb begin
        opts = _CColorizeOpts(Cstring(pointer(fb)), Cstring(pointer(cb)), lo, hi,
                              Cstring(pointer(nb)), Int64(steps),
                              component === nothing ? Int64(-1) : Int64(component),
                              Int32(symmetric), Int32(0), (0, 0, 0, 0))
        ccall(_sym(:mio_colorize), Ptr{Cvoid},
              (Ptr{Cvoid}, Ref{_CColorizeOpts}, Ref{_CColorizeReport}), _handle(m), Ref(opts),
              report)
    end
    out = Mesh(_check_ptr(ptr))
    r = report[]
    (mesh=out, vmin=r.vmin, vmax=r.vmax, location=r.location == 1 ? "cell" : "point")
end

"""
    has_netgen() -> Bool

Whether the library carries the Netgen backend [`tetrahedralize`](@ref) needs
(a build with `-DMESHIOPLUSPLUS_WITH_NETGEN=ON`).
"""
has_netgen() = ccall(_sym(:mio_has_netgen), Cint, ()) != 0

"""
    has_mmg() -> Bool

Whether the library carries the MMG backend (a build with
`-DMESHIOPLUSPLUS_WITH_MMG=ON`): `adapt`, `discretize_levelset`, `mesh_domain`
and `move` throw without it.
"""
has_mmg() = ccall(_sym(:mio_has_mmg), Cint, ()) != 0

"""
    has_mmg_elas() -> Bool

Whether the linked MMG was built with ELAS, which `move` needs.
"""
has_mmg_elas() = ccall(_sym(:mio_has_mmg_elas), Cint, ()) != 0

_cstring(s::AbstractString) = push!(Vector{UInt8}(codeunits(String(s))), 0x00)

# The shared MMG parameters, and the buffers they point into (kept alive by the caller).
function _mmg_params(mesh::Mesh, hmin, hmax, hausd, hgrad, hgradreq, angle_detection, optim,
                     noinsert, noswap, nomove, nosurf, nreg, opnbdy, mem, verbose, frozen,
                     required, local_params)
    frozen_ids = Int64[Int64(i) - 1 for i in frozen]  # 1-based in, 0-based to C
    names = [_cstring(r) for r in required]
    name_ptrs = Ptr{UInt8}[pointer(n) for n in names]
    local_names = [_cstring(String(k)) for (k, _) in local_params]
    locals = _CMmgLocal[_CMmgLocal(pointer(local_names[i]),
                                   Cdouble(get(v, "hmin", get(v, :hmin, -1.0))),
                                   Cdouble(get(v, "hmax", get(v, :hmax, -1.0))),
                                   Cdouble(get(v, "hausd", get(v, :hausd, -1.0))))
                        for (i, (_, v)) in enumerate(local_params)]
    angle = angle_detection === nothing || angle_detection === false ? -1.0 :
            Cdouble(angle_detection)
    p = _CMmgParams(Cdouble(hmin), Cdouble(hmax), Cdouble(hausd), Cdouble(hgrad),
                    Cdouble(hgradreq), angle, Int32(optim), Int32(noinsert), Int32(noswap),
                    Int32(nomove), Int32(nosurf), Int32(nreg), Int32(opnbdy), Int32(verbose),
                    Int64(mem), isempty(frozen_ids) ? Ptr{Int64}(C_NULL) : pointer(frozen_ids),
                    Int64(length(frozen_ids)),
                    isempty(name_ptrs) ? Ptr{Ptr{UInt8}}(C_NULL) : pointer(name_ptrs),
                    Int64(length(name_ptrs)),
                    isempty(locals) ? Ptr{_CMmgLocal}(C_NULL) : pointer(locals),
                    Int64(length(locals)), (0, 0, 0, 0))
    return p, (frozen_ids, names, name_ptrs, local_names, locals)
end

function _mmg_report(r::_CMmgReport)
    lib = String(UInt8[c for c in r.library if c != 0x00])
    (status=r.status == 0 ? "success" : "low_failure", library=lib,
     num_points_in=Int(r.num_points_in), num_points_out=Int(r.num_points_out),
     num_elements_in=Int(r.num_elements_in), num_elements_out=Int(r.num_elements_out),
     quality=(min=Float64(r.quality_min), mean=Float64(r.quality_mean),
              max=Float64(r.quality_max), histogram=collect(Int, r.quality_histogram)),
     num_warnings=Int(r.num_warnings))
end

"""
    adapt(mesh; hsiz=nothing, size=nothing, metric=nothing, hmin=-1, hmax=-1, hausd=-1,
          hgrad=-1, hgradreq=-1, angle_detection=45.0, optim=false, noinsert=false,
          noswap=false, nomove=false, nosurf=false, nreg=false, opnbdy=false, mem=-1,
          verbose=-1, frozen=Int[], required=String[], local_params=Dict())
        -> (mesh, status, library, num_points_in, ..., quality)

Remesh `mesh` to a size or a metric with MMG (roadmap 3.2.1, `doc/mmg.md`):
mmg3d for tetrahedra, mmg2d for a planar triangle mesh, mmgs for another
triangle surface. `size` and `metric` name point data (a metric in the
lower-triangle order of a Medit `.sol`); `frozen` lists 1-based points MMG must
keep, `required` the regions whose entities it must keep, `local_params` maps a
region to its `hmin`/`hmax`/`hausd`. Throws naming `MESHIOPLUSPLUS_WITH_MMG`
without the backend ([`has_mmg`](@ref)).
"""
function adapt(mesh::Mesh; hsiz=nothing, size=nothing, metric=nothing, hmin::Real=-1.0,
               hmax::Real=-1.0, hausd::Real=-1.0, hgrad::Real=-1.0, hgradreq::Real=-1.0,
               angle_detection=45.0, optim::Bool=false, noinsert::Bool=false,
               noswap::Bool=false, nomove::Bool=false, nosurf::Bool=false, nreg::Bool=false,
               opnbdy::Bool=false, mem::Integer=-1, verbose::Integer=-1,
               frozen=Int[], required=String[], local_params=Dict{String,Any}())
    params, keep = _mmg_params(mesh, hmin, hmax, hausd, hgrad, hgradreq, angle_detection,
                               optim, noinsert, noswap, nomove, nosurf, nreg, opnbdy, mem,
                               verbose, frozen, required, local_params)
    size_c = size === nothing ? UInt8[0x00] : _cstring(size)
    metric_c = metric === nothing ? UInt8[0x00] : _cstring(metric)
    report = Ref{_CMmgReport}()
    ptr = GC.@preserve keep size_c metric_c begin
        opts = _CAdaptOpts(params, hsiz === nothing ? -1.0 : Cdouble(hsiz),
                           size === nothing ? Ptr{UInt8}(C_NULL) : pointer(size_c),
                           metric === nothing ? Ptr{UInt8}(C_NULL) : pointer(metric_c),
                           (0, 0, 0, 0))
        ccall(_sym(:mio_adapt), Ptr{Cvoid}, (Ptr{Cvoid}, Ref{_CAdaptOpts}, Ref{_CMmgReport}),
              _handle(mesh), Ref(opts), report)
    end
    out = Mesh(_check_ptr(ptr))
    Base.merge((mesh=out,), _mmg_report(report[]))
end

"""
    discretize_levelset(mesh, level_set; isovalue=0.0, remove_small_components=nothing,
                        materials=[], <adapt's MMG parameters>)
        -> (mesh, status, library, num_points_in, ..., quality)

Embed the `isovalue` of the one-component point data `level_set` in `mesh`
conformally and remesh both sides with MMG (roadmap 3.2.2, `doc/mmg.md`). The
elements carry cell data `medit:ref` (3 inside, 2 outside) and the cell regions
`inside`/`outside`; the interface is the side region `interface`.
`remove_small_components` is MMG's `rmc` volume fraction; `materials` lists
`(ref, split, ref_inside, ref_outside)` tuples for MMG's multi-material mode.
"""
function discretize_levelset(mesh::Mesh, level_set::AbstractString; isovalue::Real=0.0,
                             remove_small_components=nothing, materials=(),
                             hmin::Real=-1.0, hmax::Real=-1.0, hausd::Real=-1.0,
                             hgrad::Real=-1.0, hgradreq::Real=-1.0, angle_detection=45.0,
                             optim::Bool=false, noinsert::Bool=false, noswap::Bool=false,
                             nomove::Bool=false, nosurf::Bool=false, nreg::Bool=false,
                             opnbdy::Bool=false, mem::Integer=-1, verbose::Integer=-1,
                             frozen=Int[], required=String[], local_params=Dict{String,Any}())
    params, keep = _mmg_params(mesh, hmin, hmax, hausd, hgrad, hgradreq, angle_detection,
                               optim, noinsert, noswap, nomove, nosurf, nreg, opnbdy, mem,
                               verbose, frozen, required, local_params)
    ls_c = _cstring(level_set)
    mats = _CMmgMaterial[_CMmgMaterial(Int64(m[1]), Int32(Bool(m[2])), Int32(0),
                                       Int64(m[3]), Int64(m[4])) for m in materials]
    report = Ref{_CMmgReport}()
    ptr = GC.@preserve keep ls_c mats begin
        opts = _CDiscretizeLevelsetOpts(params, pointer(ls_c), Cdouble(isovalue),
                                        remove_small_components === nothing ? -1.0 :
                                        Cdouble(remove_small_components),
                                        isempty(mats) ? Ptr{_CMmgMaterial}(C_NULL) :
                                        pointer(mats), Int64(length(mats)), (0, 0, 0, 0))
        ccall(_sym(:mio_discretize_levelset), Ptr{Cvoid},
              (Ptr{Cvoid}, Ref{_CDiscretizeLevelsetOpts}, Ref{_CMmgReport}),
              _handle(mesh), Ref(opts), report)
    end
    out = Mesh(_check_ptr(ptr))
    Base.merge((mesh=out,), _mmg_report(report[]))
end

"""
    mesh_domain(boundary; hsiz=nothing, size=nothing, holes=[], <adapt's MMG parameters>)
        -> (mesh, num_subdomains, num_holes, status, library, ...)

Triangulate the planar domain bounded by the `line` cells of `boundary` with
MMG (`MMG2D_mmg2dmesh`, roadmap 3.2.3, `doc/mmg.md`). Each enclosed region is a
cell region `subdomain:<k>`; `holes` lists one `(x, y)` point inside each region
to leave out.
"""
function mesh_domain(boundary::Mesh; hsiz=nothing, size=nothing, holes=(),
                     hmin::Real=-1.0, hmax::Real=-1.0, hausd::Real=-1.0, hgrad::Real=-1.0,
                     hgradreq::Real=-1.0, angle_detection=45.0, optim::Bool=false,
                     noinsert::Bool=false, noswap::Bool=false, nomove::Bool=false,
                     nosurf::Bool=false, nreg::Bool=false, opnbdy::Bool=false,
                     mem::Integer=-1, verbose::Integer=-1, frozen=Int[], required=String[],
                     local_params=Dict{String,Any}())
    params, keep = _mmg_params(boundary, hmin, hmax, hausd, hgrad, hgradreq, angle_detection,
                               optim, noinsert, noswap, nomove, nosurf, nreg, opnbdy, mem,
                               verbose, frozen, required, local_params)
    size_c = size === nothing ? UInt8[0x00] : _cstring(size)
    xy = Cdouble[c for h in holes for c in (Cdouble(h[1]), Cdouble(h[2]))]
    report = Ref{_CMmgReport}()
    nsub = Ref{Int64}(0)
    nholes = Ref{Int64}(0)
    ptr = GC.@preserve keep size_c xy begin
        opts = _CMeshDomainOpts(params, hsiz === nothing ? -1.0 : Cdouble(hsiz),
                                size === nothing ? Ptr{UInt8}(C_NULL) : pointer(size_c),
                                isempty(xy) ? Ptr{Cdouble}(C_NULL) : pointer(xy),
                                Int64(length(xy) ÷ 2), (0, 0, 0, 0))
        ccall(_sym(:mio_mesh_domain), Ptr{Cvoid},
              (Ptr{Cvoid}, Ref{_CMeshDomainOpts}, Ref{_CMmgReport}, Ref{Int64}, Ref{Int64}),
              _handle(boundary), Ref(opts), report, nsub, nholes)
    end
    out = Mesh(_check_ptr(ptr))
    Base.merge((mesh=out, num_subdomains=Int(nsub[]), num_holes=Int(nholes[])),
          _mmg_report(report[]))
end

"""
    move(mesh, displacement, moving; lagrangian=1, <adapt's MMG parameters>)
        -> (mesh, status, library, num_points_in, ..., quality)

Move `mesh` by the point-data vector `displacement` of the boundary facets of
the regions `moving` (a name or a list) and remesh as it deforms, with MMG's
Lagrangian mode (roadmap 3.2.4, `doc/mmg.md`): MMG extends the displacement
inward by linear elasticity (ELAS) and keeps the rest of the boundary in place.
`lagrangian` is MMG's `lag` (0, 1 or 2). Throws naming `MESHIOPLUSPLUS_WITH_MMG`
without MMG, or ELAS without it ([`has_mmg_elas`](@ref)).
"""
function move(mesh::Mesh, displacement::AbstractString, moving; lagrangian::Integer=1,
              hmin::Real=-1.0, hmax::Real=-1.0, hausd::Real=-1.0, hgrad::Real=-1.0,
              hgradreq::Real=-1.0, angle_detection=45.0, optim::Bool=false,
              noinsert::Bool=false, noswap::Bool=false, nomove::Bool=false, nosurf::Bool=false,
              nreg::Bool=false, opnbdy::Bool=false, mem::Integer=-1, verbose::Integer=-1,
              frozen=Int[], required=String[], local_params=Dict{String,Any}())
    params, keep = _mmg_params(mesh, hmin, hmax, hausd, hgrad, hgradreq, angle_detection,
                               optim, noinsert, noswap, nomove, nosurf, nreg, opnbdy, mem,
                               verbose, frozen, required, local_params)
    disp_c = _cstring(displacement)
    names = [_cstring(r) for r in (moving isa AbstractString ? [moving] : moving)]
    name_ptrs = Ptr{UInt8}[pointer(n) for n in names]
    report = Ref{_CMmgReport}()
    ptr = GC.@preserve keep disp_c names name_ptrs begin
        opts = _CMoveOpts(params, pointer(disp_c),
                          isempty(name_ptrs) ? Ptr{Ptr{UInt8}}(C_NULL) : pointer(name_ptrs),
                          Int64(length(name_ptrs)), Int32(lagrangian), Int32(0), (0, 0, 0, 0))
        ccall(_sym(:mio_move), Ptr{Cvoid}, (Ptr{Cvoid}, Ref{_CMoveOpts}, Ref{_CMmgReport}),
              _handle(mesh), Ref(opts), report)
    end
    out = Mesh(_check_ptr(ptr))
    Base.merge((mesh=out,), _mmg_report(report[]))
end

"""
    metric_from_hessian(mesh, hessian; error=1e-2, hmin=1e-3, hmax=1.0, output="metric") -> Mesh

The anisotropic metric of the row-major point-data Hessian `hessian` ((n, 9), or
(n, 4) in 2-D) for [`adapt`](@ref): `M = c/error |H|`, its sizes clamped to
`[hmin, hmax]`, as point data `output` in the lower-triangle order of a Medit
`.sol`. Needs no MMG.
"""
function metric_from_hessian(mesh::Mesh, hessian::AbstractString; error::Real=1e-2,
                             hmin::Real=1e-3, hmax::Real=1.0, output::AbstractString="metric")
    Mesh(_check_ptr(ccall(_sym(:mio_metric_from_hessian), Ptr{Cvoid},
                          (Ptr{Cvoid}, Cstring, Cdouble, Cdouble, Cdouble, Cstring),
                          _handle(mesh), hessian, error, hmin, hmax, output)))
end

"""
    tetrahedralize(mesh; backend="auto", max_size=0.0, max_volume=0.0, grading=0.3,
                   optimize_steps=3, max_cells=50_000_000, face_array="cad:face")
        -> (; mesh, num_tetrahedra, num_points, num_new_points, num_regions,
              num_boundary_triangles, boundary_preserved, hausdorff_to_surface,
              min_dihedral_deg, max_dihedral_deg, volume, backend)

Fill a closed, manifold, consistently oriented triangle surface with tetrahedra
(roadmap 6.7): through Netgen its triangles and points are kept exactly, through
the native octree backend its boundary is fitted within `hausdorff_to_surface`.
Each shell bounds a region (a box inside a box gives two), meshed on its own.
`max_size` and `max_volume` (0: none) bound the tetrahedra, `grading` is how fast
the size grows inward, `optimize_steps` the backend's volume optimisation
passes; `face_array` names the CAD face cell data carried to the triangles.
The result has a `tetra` block, a `triangle` block (the surface's with Netgen,
its input points first; the tetrahedra's boundary with the octree backend) and
the cell data `tetrahedralize:region`. `backend` is `"auto"` (Netgen when
[`has_netgen`](@ref), else `"octree"`), `"netgen"` (throws naming
`MESHIOPLUSPLUS_WITH_NETGEN` without it), `"octree"`, or another name of
[`tetrahedralize_backends`](@ref): `"tetgen"` and `"gmsh"` run those installed
programs as separate processes, the surface kept (checked).
See `doc/tetrahedralize.md`.
"""
function tetrahedralize(mesh::Mesh; backend::AbstractString="auto", max_size::Real=0.0,
                        max_volume::Real=0.0, grading::Real=0.3, optimize_steps::Integer=3,
                        max_cells::Integer=50_000_000, face_array::AbstractString="cad:face")
    isempty(backend) &&
        throw(ArgumentError("meshio++: tetrahedralize: backend must be \"auto\", \"netgen\", \"octree\" or a mesher's name"))
    # Any other name is a named mesher (3); the library refuses an unknown one,
    # listing the available names.
    backend_c = backend == "auto" ? Int32(0) : backend == "netgen" ? Int32(1) :
                backend == "octree" ? Int32(2) : Int32(3)
    fname = Vector{UInt8}(codeunits(String(face_array)))
    push!(fname, 0x00)
    pname = Vector{UInt8}(codeunits(String(backend)))
    push!(pname, 0x00)
    report = Ref{_CTetrahedralizeReport}()
    ptr = GC.@preserve fname pname begin
        opts = _CTetrahedralizeOpts(Cdouble(max_size), Cdouble(max_volume), Cdouble(grading),
                                    Int64(max_cells), pointer(fname), backend_c,
                                    Int32(optimize_steps),
                                    backend_c == 3 ? pointer(pname) : Ptr{UInt8}(C_NULL),
                                    (0, 0, 0))
        ccall(_sym(:mio_tetrahedralize), Ptr{Cvoid},
              (Ptr{Cvoid}, Ref{_CTetrahedralizeOpts}, Ref{_CTetrahedralizeReport}),
              _handle(mesh), Ref(opts), report)
    end
    out = Mesh(_check_ptr(ptr))
    r = report[]
    (mesh=out, num_tetrahedra=Int(r.num_tetrahedra), num_points=Int(r.num_points),
     num_new_points=Int(r.num_new_points), num_regions=Int(r.num_regions),
     num_boundary_triangles=Int(r.num_boundary_triangles),
     boundary_preserved=r.boundary_preserved != 0,
     hausdorff_to_surface=Float64(r.hausdorff_to_surface),
     min_dihedral_deg=Float64(r.min_dihedral_deg), max_dihedral_deg=Float64(r.max_dihedral_deg),
     volume=Float64(r.volume),
     backend=r.backend_name == C_NULL ? (r.backend == 1 ? "netgen" : "octree") :
             unsafe_string(r.backend_name))
end

"""
    remesh_quads(mesh; target_faces=0, backend="auto", sharp=false, boundary=true,
                 adaptive=false, dominant=false, seed=0, max_cells=50_000_000)
        -> (; mesh, num_quads, num_triangles, num_points, min_scaled_jacobian,
              mean_scaled_jacobian, hausdorff_to_surface, backend)

Remesh a triangle surface into field-aligned quads (roadmap 7.7.3) through an
installed QuadriFlow or Instant Meshes, run as a separate process. `target_faces`
0 takes half the triangles; `backend` is `"auto"` (QuadriFlow, else Instant
Meshes) or a name of [`remesh_quads_backends`](@ref). No remesher installed
throws a `MeshioError` naming the programs. See `doc/remesh_quads.md`.
"""
function remesh_quads(mesh::Mesh; target_faces::Integer=0, backend::AbstractString="auto",
                      sharp::Bool=false, boundary::Bool=true, adaptive::Bool=false,
                      dominant::Bool=false, seed::Integer=0, max_cells::Integer=50_000_000)
    bname = Vector{UInt8}(codeunits(String(backend)))
    push!(bname, 0x00)
    report = Ref{_CRemeshQuadsReport}()
    ptr = GC.@preserve bname begin
        opts = _CRemeshQuadsOpts(pointer(bname), Int64(target_faces), Int64(max_cells),
                                 Int32(sharp), Int32(boundary), Int32(adaptive),
                                 Int32(dominant), Int32(seed), Int32(0), (0, 0, 0, 0))
        ccall(_sym(:mio_remesh_quads), Ptr{Cvoid},
              (Ptr{Cvoid}, Ref{_CRemeshQuadsOpts}, Ref{_CRemeshQuadsReport}),
              _handle(mesh), Ref(opts), report)
    end
    out = Mesh(_check_ptr(ptr))
    r = report[]
    (mesh=out, num_quads=Int(r.num_quads), num_triangles=Int(r.num_triangles),
     num_points=Int(r.num_points), min_scaled_jacobian=Float64(r.min_scaled_jacobian),
     mean_scaled_jacobian=Float64(r.mean_scaled_jacobian),
     hausdorff_to_surface=Float64(r.hausdorff_to_surface),
     backend=r.backend_name == C_NULL ? "" : unsafe_string(r.backend_name))
end

"""
    has_cutcells() -> Bool

Whether the library carries CutCells, which [`cut_level_set`](@ref) and
[`cut_quadrature`](@ref) need (a build with `MESHIOPLUSPLUS_WITH_CUTCELLS`).
"""
has_cutcells() = ccall(_sym(:mio_has_cutcells), Int32, ()) != 0

function _cut_call(sym::Symbol, mesh::Mesh, level_sets, planes, selection, include_uncut,
                   backend, degree, triangulate, order, record_parent_ids)
    names = [Vector{UInt8}(codeunits(String(a) * "\0")) for a in level_sets]
    ptrs = [Cstring(pointer(n)) for n in names]
    flat = Cdouble[]
    for (origin, normal) in planes
        length(origin) == 3 && length(normal) == 3 ||
            throw(ArgumentError("meshio++: cut: a plane is (origin, normal), three numbers each"))
        append!(flat, Cdouble.(origin))
        append!(flat, Cdouble.(normal))
    end
    sel = Vector{UInt8}(codeunits(String(selection) * "\0"))
    bk = Vector{UInt8}(codeunits(String(backend) * "\0"))
    report = Ref{_CCutReport}()
    ptr = GC.@preserve names ptrs flat sel bk begin
        opts = _CCutOpts(isempty(ptrs) ? Ptr{Cstring}(C_NULL) : pointer(ptrs), Int64(length(ptrs)),
                         isempty(flat) ? Ptr{Cdouble}(C_NULL) : pointer(flat),
                         Int64(length(flat) ÷ 6), pointer(sel), pointer(bk), Int32(include_uncut),
                         Int32(degree), Int32(triangulate), Int32(order), Int32(record_parent_ids),
                         Int32(0), (0, 0, 0, 0))
        ccall(_sym(sym), Ptr{Cvoid}, (Ptr{Cvoid}, Ref{_CCutOpts}, Ref{_CCutReport}),
              _handle(mesh), Ref(opts), report)
    end
    Mesh(_check_ptr(ptr)), report[]
end

"""
    cut_level_set(mesh; level_sets=String[], planes=[], selection="", include_uncut=true,
                  backend="auto", degree=2, triangulate=false, record_parent_ids=true)
        -> (; mesh, num_cells, num_cut_cells, num_uncut_cells, num_pieces, dimension,
              measure, backend)

The part of `mesh` that level sets select, through CutCells (`doc/cut.md`):
`level_sets` are scalar point-data arrays, `planes` `(origin, normal)` pairs
named `plane0`, `plane1`, ...; `selection` an expression such as `"phi < 0"`,
`"phi = 0"` (the interface) or `"a < 0 and b > 0"` (empty: `"<first> < 0"`).
`backend` is `"lut"` (straight pieces), `"quadrays"` (curved pieces of
`degree`) or `"auto"`. Without CutCells it throws a `MeshioError` naming
`MESHIOPLUSPLUS_WITH_CUTCELLS`.
"""
function cut_level_set(mesh::Mesh; level_sets=String[], planes=(), selection::AbstractString="",
                       include_uncut::Bool=true, backend::AbstractString="auto",
                       degree::Integer=2, triangulate::Bool=false, record_parent_ids::Bool=true)
    out, r = _cut_call(:mio_cut_level_set, mesh, level_sets, planes, selection, include_uncut,
                       backend, degree, triangulate, 2, record_parent_ids)
    (mesh=out, num_cells=Int(r.num_cells), num_cut_cells=Int(r.num_cut_cells),
     num_uncut_cells=Int(r.num_uncut_cells), num_pieces=Int(r.num_pieces),
     dimension=Int(r.dimension), measure=Float64(r.measure),
     backend=r.backend_name == C_NULL ? "" : unsafe_string(r.backend_name))
end

"""
    cut_quadrature(mesh; level_sets=String[], planes=[], selection="", order=2, ...)
        -> (; mesh, num_rules, num_points, sum_weights, dimension, backend)

Quadrature rules of degree `order` on the part [`cut_level_set`](@ref) would
give: a point cloud with point data `quadrature:weight`, `quadrature:reference`,
`quadrature:rule` and `quadrature:parent_cell`.
"""
function cut_quadrature(mesh::Mesh; level_sets=String[], planes=(), selection::AbstractString="",
                        order::Integer=2, include_uncut::Bool=true,
                        backend::AbstractString="auto", triangulate::Bool=false,
                        record_parent_ids::Bool=true)
    out, r = _cut_call(:mio_cut_quadrature, mesh, level_sets, planes, selection, include_uncut,
                       backend, 2, triangulate, order, record_parent_ids)
    (mesh=out, num_rules=Int(r.num_rules), num_points=Int(r.num_points),
     sum_weights=Float64(r.measure), dimension=Int(r.dimension),
     backend=r.backend_name == C_NULL ? "" : unsafe_string(r.backend_name))
end

"""
    remesh_quads_backends() -> Vector{String}

The quad remeshers [`remesh_quads`](@ref) takes by name besides `"auto"`, sorted:
the installed `"quadriflow"` and `"instant_meshes"` programs, found on `PATH`.
"""
function remesh_quads_backends()
    n = ccall(_sym(:mio_quad_remesher_count), Int64, ())
    n < 0 && throw(MeshioError(MIO_ERR_INTERNAL, last_error()))
    [unsafe_string(ccall(_sym(:mio_quad_remesher_name), Ptr{UInt8}, (Int64,), i))
     for i in 0:(n - 1)]
end

"""
    tetrahedralize_backends() -> Vector{String}

The `backend` names [`tetrahedralize`](@ref) takes besides `"auto"`, sorted:
`"octree"`, then every registered mesher (`"netgen"` when built in, the installed
`"tetgen"` and `"gmsh"`, found on `PATH`).
"""
function tetrahedralize_backends()
    n = ccall(_sym(:mio_tet_mesher_count), Int64, ())
    n < 0 && throw(MeshioError(MIO_ERR_INTERNAL, last_error()))
    [unsafe_string(ccall(_sym(:mio_tet_mesher_name), Ptr{UInt8}, (Int64,), i))
     for i in 0:(n - 1)]
end

"""
    interpolate(source, target; method="nearest", arrays=String[],
                extrapolate=false, default_value=0.0, on_conflict="error") -> Mesh

Sample data arrays from `source` onto `target` (cross-mesh field transfer).
The result is a copy of `target` — geometry, connectivity and its own data
preserved exactly — with the requested source arrays sampled onto it.

An empty `arrays` means every source point-data array; cell data transfers only
when named explicitly. `method` is `"nearest"` or `"barycentric"`;
`on_conflict` is `"error"`, `"overwrite"` or `"suffix"`. See
`doc/interpolate.md`.
"""
function interpolate(source::Mesh, target::Mesh; method::AbstractString="nearest",
                     arrays=String[], extrapolate::Bool=false,
                     default_value::Real=0.0, on_conflict::AbstractString="error")
    sh = _handle(source); th = _handle(target)
    ptr = _with_names(arrays) do names_ptr, count
        ccall(_sym(:mio_interpolate), Ptr{Cvoid},
              (Ptr{Cvoid}, Ptr{Cvoid}, Cstring, Ptr{Cstring}, Int64, Cint, Cdouble, Cstring),
              sh, th, method, names_ptr, count, extrapolate ? Cint(1) : Cint(0),
              Float64(default_value), on_conflict)
    end
    Mesh(_check_ptr(ptr))
end

"""
    conservative_interpolate(source, target; arrays=String[],
                             default_value=0.0, on_conflict="error") -> Mesh

Mass-preserving cross-mesh field transfer: an exact overlap-measure weighted
remap, so that over the region the two meshes share, `sum(target value *
target measure)` equals `sum(source value * source measure)` — the property
[`interpolate`](@ref)'s `"barycentric"` mode does not have. Both meshes are
simplexified first (accepting ragged/polyhedron blocks for free).

Unlike `interpolate`, an empty `arrays` means every source point_data AND
cell_data array — there is one algorithm regardless of location. Output
arrays are always Float64. `on_conflict` is `"error"`, `"overwrite"` or
`"suffix"`. See `doc/conservative_interpolate.md`.
"""
function conservative_interpolate(source::Mesh, target::Mesh; arrays=String[],
                                   default_value::Real=0.0,
                                   on_conflict::AbstractString="error")
    sh = _handle(source); th = _handle(target)
    ptr = _with_names(arrays) do names_ptr, count
        ccall(_sym(:mio_conservative_interpolate), Ptr{Cvoid},
              (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cstring}, Int64, Cdouble, Cstring),
              sh, th, names_ptr, count, Float64(default_value), on_conflict)
    end
    Mesh(_check_ptr(ptr))
end

"""
    meshes_equal(a, b; atol=0.0, rtol=0.0, unordered=false) -> Bool

Are two meshes equal within tolerance (i.e. the diff verdict is not
`:different`)? Named regions are **not** considered — they never reach the C++
comparison.
"""
function meshes_equal(a::Mesh, b::Mesh; atol::Real=0.0, rtol::Real=0.0,
                      unordered::Bool=false)
    out = Ref{Cint}(0)
    _check(ccall(_sym(:mio_meshes_equal), Cint,
                 (Ptr{Cvoid}, Ptr{Cvoid}, Cdouble, Cdouble, Cint, Ptr{Cint}),
                 _handle(a), _handle(b), Float64(atol), Float64(rtol),
                 unordered ? Cint(1) : Cint(0), out))
    out[] != 0
end

"""
    DiffReport

Structured result of [`diff`](@ref). `verdict` is `:identical`,
`:equal_within_tolerance` or `:different`.
"""
struct DiffReport
    verdict::Symbol
    points::NamedTuple{(:max_abs_error, :max_rel_error, :worst_index, :num_exceeding),
                       Tuple{Float64,Float64,Int,Int}}
    blocks::Vector{NamedTuple{(:type_mismatch, :count_mismatch, :conn_mismatch_count),
                              Tuple{Bool,Bool,Int}}}
end

const _VERDICTS = (:identical, :equal_within_tolerance, :different)

"""
    diff(a, b; atol=0.0, rtol=0.0, unordered=false) -> DiffReport

Compare two meshes: points, cell blocks and point/cell/field data, within
`abs_err <= atol + rtol*|expected|` with `expected` taken from `a`.
`unordered` matches points by spatial proximity instead of by index.

Named regions are **not** compared (they never reach the C++ core). Not
exported (it would shadow `Base.diff`): call it qualified.
"""
function diff(a::Mesh, b::Mesh; atol::Real=0.0, rtol::Real=0.0, unordered::Bool=false)
    out = Ref{Ptr{Cvoid}}(C_NULL)
    _check(ccall(_sym(:mio_diff), Cint,
                 (Ptr{Cvoid}, Ptr{Cvoid}, Cdouble, Cdouble, Cint, Ptr{Ptr{Cvoid}}),
                 _handle(a), _handle(b), Float64(atol), Float64(rtol),
                 unordered ? Cint(1) : Cint(0), out))
    result = out[]
    try
        v = ccall(_sym(:mio_diff_result_verdict), Cint, (Ptr{Cvoid},), result)
        mabs = Ref{Cdouble}(0.0); mrel = Ref{Cdouble}(0.0)
        worst = Ref{Int64}(0); nex = Ref{Int64}(0)
        _check(ccall(_sym(:mio_diff_result_point_summary), Cint,
                     (Ptr{Cvoid}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Int64}, Ptr{Int64}),
                     result, mabs, mrel, worst, nex))
        nb = ccall(_sym(:mio_diff_result_num_block_diffs), Int64, (Ptr{Cvoid},), result)
        nb < 0 && (nb = 0)
        blocks = NamedTuple{(:type_mismatch, :count_mismatch, :conn_mismatch_count),
                            Tuple{Bool,Bool,Int}}[]
        for i in 0:(nb-1)
            tm = Ref{Cint}(0); cm = Ref{Cint}(0); cc = Ref{Int64}(0)
            _check(ccall(_sym(:mio_diff_result_block), Cint,
                         (Ptr{Cvoid}, Int64, Ptr{Cint}, Ptr{Cint}, Ptr{Int64}),
                         result, i, tm, cm, cc))
            push!(blocks, (type_mismatch=tm[] != 0, count_mismatch=cm[] != 0,
                           conn_mismatch_count=Int(cc[])))
        end
        # worst_index is a flat index into the coordinate buffer: 1-based here,
        # with the C API's -1 "no offender" becoming 0.
        w = worst[] < 0 ? 0 : Int(worst[]) + 1
        DiffReport(_VERDICTS[v+1],
                   (max_abs_error=Float64(mabs[]), max_rel_error=Float64(mrel[]),
                    worst_index=w, num_exceeding=Int(nex[])), blocks)
    finally
        ccall(_sym(:mio_diff_result_free), Cvoid, (Ptr{Cvoid},), result)
    end
end

# --- statistics --------------------------------------------------------------

"""
    stats(mesh) -> NamedTuple

Read-only geometric statistics: counts, bounding box, extent, centroid, total
area, signed / unsigned volume and the inverted-cell count.

Per-cell-type counts are **not carried across the C ABI** (a documented gap);
use [`cell_block_types`](@ref) with [`cell_block_info`](@ref) instead.
"""
function stats(m::Mesh)
    out = Ref{_CStatsReport}()
    _check(ccall(_sym(:mio_stats), Cint, (Ptr{Cvoid}, Ptr{_CStatsReport}), _handle(m), out))
    r = out[]
    (num_points=Int(r.num_points), num_cells=Int(r.num_cells),
     bbox_min=r.bbox_min, bbox_max=r.bbox_max, extent=r.extent, centroid=r.centroid,
     total_area=r.total_area, signed_volume=r.signed_volume,
     unsigned_volume=r.unsigned_volume, num_inverted=Int(r.num_inverted))
end

"""
    compute_bandwidth(mesh) -> Int

Connectivity bandwidth: the maximum over all cells of (max node index − min
node index) within a cell. A before/after measure for [`reorder`](@ref).
"""
compute_bandwidth(m::Mesh) =
    Int(_check_count(ccall(_sym(:mio_compute_bandwidth), Int64, (Ptr{Cvoid},), _handle(m)),
                     "mio_compute_bandwidth"))

# --- opaque-result operations ------------------------------------------------

"""
    reorder(mesh, method) -> (; mesh, node_perm, cell_perms)

Renumber a mesh as a **pure permutation** of node and element indices, to
reduce sparse-matrix bandwidth or improve cache locality. `method` is `"rcm"`,
`"morton"` or `"hilbert"`.

`node_perm` and each entry of `cell_perms` are 1-based old→new index vectors.
See `doc/reorder.md`.
"""
function reorder(m::Mesh, method::AbstractString)
    result = _check_ptr(ccall(_sym(:mio_reorder), Ptr{Cvoid}, (Ptr{Cvoid}, Cstring),
                              _handle(m), method))
    try
        node_perm = _result_map(result, :mio_reorder_result_node_perm)
        cell_perms = _result_cell_maps(result, :mio_reorder_result_num_cell_perms,
                                       :mio_reorder_result_cell_perm)
        out = Mesh(_check_ptr(ccall(_sym(:mio_reorder_result_take_mesh), Ptr{Cvoid},
                                    (Ptr{Cvoid},), result)))
        (mesh=out, node_perm=node_perm, cell_perms=cell_perms)
    finally
        ccall(_sym(:mio_reorder_result_free), Cvoid, (Ptr{Cvoid},), result)
    end
end

"""
    split(mesh; by="type", tag="") -> Vector{Pair{String,Mesh}}

Partition a mesh into submeshes by a criterion: `by` is `"type"` (one piece per
cell type), `"component"` (connected components) or `"region"`/`"tag"` (the
distinct values of an integer cell-data array, named by `tag`, or the first
integer array when `tag` is empty).

Each pair is `key => piece`. Every piece owns its handle. Not exported (it
would shadow `Base.split`): call it qualified. See `doc/split.md`.
"""
function split(m::Mesh; by::AbstractString="type", tag::AbstractString="")
    result = _check_ptr(ccall(_sym(:mio_split), Ptr{Cvoid}, (Ptr{Cvoid}, Cstring, Cstring),
                              _handle(m), by, tag))
    try
        n = _check_count(ccall(_sym(:mio_split_result_count), Int64, (Ptr{Cvoid},), result),
                         "mio_split_result_count")
        out = Pair{String,Mesh}[]
        for i in 0:(n-1)
            key = _getstring() do buf, len
                ccall(_sym(:mio_split_result_key), Int64,
                      (Ptr{Cvoid}, Int64, Ptr{UInt8}, Int64), result, i, buf, len)
            end
            piece = Mesh(_check_ptr(ccall(_sym(:mio_split_result_take_mesh), Ptr{Cvoid},
                                          (Ptr{Cvoid}, Int64), result, i)))
            push!(out, key => piece)
        end
        out
    finally
        ccall(_sym(:mio_split_result_free), Cvoid, (Ptr{Cvoid},), result)
    end
end

"""
    convert_cells(mesh, mode; record_parent_ids=false) -> (; mesh, point_map, cell_maps)

Convert the **element representation** (not the file format). `mode` is
`"linearize"` (higher-order cells → their linear base, pruning orphaned
nodes), `"simplexify"` (decompose into same-dimension simplices) or
`"elevate"` (linear → serendipity quadratic).

`point_map` is 1-based input-point → output-point, `0` where the point was
pruned; each `cell_maps[b]` is input cell → its **first child** in the
corresponding output block. See `doc/convert_cells.md`.
"""
function convert_cells(m::Mesh, mode::AbstractString; record_parent_ids::Bool=false)
    result = _check_ptr(ccall(_sym(:mio_convert_cells), Ptr{Cvoid},
                              (Ptr{Cvoid}, Cstring, Cint), _handle(m), mode,
                              record_parent_ids ? Cint(1) : Cint(0)))
    try
        pm = _result_map(result, :mio_convert_cells_result_point_map)
        cm = _result_cell_maps(result, :mio_convert_cells_result_num_cell_maps,
                               :mio_convert_cells_result_cell_map)
        out = Mesh(_check_ptr(ccall(_sym(:mio_convert_cells_result_take_mesh), Ptr{Cvoid},
                                    (Ptr{Cvoid},), result)))
        (mesh=out, point_map=pm, cell_maps=cm)
    finally
        ccall(_sym(:mio_convert_cells_result_free), Cvoid, (Ptr{Cvoid},), result)
    end
end

"""
    subdivide(mesh; record_parent_ids=false) -> (; mesh, cell_maps)

Polyhedrally refine: split every eligible 3-D cell into one polyhedral child
per face, connected to a new interior point. Needs no per-type template
table -- tabulated types (reduced to corners for a quadratic variant) and
existing polyhedron blocks are handled uniformly. Automatically conforming,
unlike [`refine`](@ref): a shared face between two input cells is never
touched. Non-3-D blocks and the full-Lagrange family (no face table) pass
through unchanged.

Unlike [`convert_cells`](@ref), there is no `point_map` -- subdivide never
prunes or renumbers an original point. Each `cell_maps[b]` is input cell →
the index of its **first** child (one per face) in the corresponding output
block. See `doc/subdivide.md`.
"""
function subdivide(m::Mesh; record_parent_ids::Bool=false)
    result = _check_ptr(ccall(_sym(:mio_subdivide), Ptr{Cvoid},
                              (Ptr{Cvoid}, Cint), _handle(m),
                              record_parent_ids ? Cint(1) : Cint(0)))
    try
        cm = _result_cell_maps(result, :mio_subdivide_result_num_cell_maps,
                               :mio_subdivide_result_cell_map)
        out = Mesh(_check_ptr(ccall(_sym(:mio_subdivide_result_take_mesh), Ptr{Cvoid},
                                    (Ptr{Cvoid},), result)))
        (mesh=out, cell_maps=cm)
    finally
        ccall(_sym(:mio_subdivide_result_free), Cvoid, (Ptr{Cvoid},), result)
    end
end

"""
    agglomerate(mesh; target_group_size=8) -> (; mesh, cell_map)

Polyhedrally coarsen: merge groups of cells into single larger polyhedral
cells via greedy seed-and-grow over the mesh's shared-face dual, absorbing
face-adjacent neighbours into a group until it reaches `target_group_size`
(a short group at a mesh boundary or pocket is expected, not an error).
Non-volume blocks pass through unchanged; points are never pruned or
renumbered ([`clean`](@ref) with `remove_orphans=true` is the follow-up for
a minimal point set). `target_group_size=1` groups every cell by itself.

Unlike [`subdivide`](@ref)'s per-block `cell_maps`, `cell_map` is a single
**flat** 1-based array (input global cell → output global cell) -- an
output cell's index is a function of which group it joined, not which
input block it came from. See `doc/agglomerate.md`.
"""
function agglomerate(m::Mesh; target_group_size::Integer=8, merge_coplanar_faces::Bool=false,
                     coplanar_angle::Real=1.0, min_sphericity::Real=0.0)
    opts = _CAgglomerateOpts(Int64(target_group_size), Int32(merge_coplanar_faces), Int32(0),
                             Float64(coplanar_angle), Float64(min_sphericity),
                             ntuple(_ -> Int64(0), 6))
    merged, rejected = Ref{Int64}(0), Ref{Int64}(0)
    result = _check_ptr(ccall(_sym(:mio_agglomerate_ex), Ptr{Cvoid},
                              (Ptr{Cvoid}, Ref{_CAgglomerateOpts}, Ptr{Int64}, Ptr{Int64}),
                              _handle(m), Ref(opts), merged, rejected))
    try
        cm = _result_map(result, :mio_agglomerate_result_cell_map)
        out = Mesh(_check_ptr(ccall(_sym(:mio_agglomerate_result_take_mesh), Ptr{Cvoid},
                                    (Ptr{Cvoid},), result)))
        (mesh=out, cell_map=cm, num_faces_merged=Int(merged[]), num_rejected=Int(rejected[]))
    finally
        ccall(_sym(:mio_agglomerate_result_free), Cvoid, (Ptr{Cvoid},), result)
    end
end

"""
    tet_to_poly(mesh; feature_angle=1.0, point_data_to_cell_data=false,
                record_vertex_ids=false) -> (; mesh, vertex_map, num_new_points,
                                              num_boundary_faces, num_split_vertices)

Convert a tetrahedral mesh into its polyhedral (median) dual: one polyhedron
per vertex, bounded by the ring of tetrahedron and triangle centroids around
each edge leaving it and, on the boundary, by its share of the boundary
triangles, fused while their normals differ by at most `feature_angle`
degrees (so feature edges and corners survive). A vertex on an interface
between tetra blocks gets one cell per block. Each tetra block becomes a
polyhedron block and each triangle block a polygon block; any other block is
dropped. The input points keep their indices.

`vertex_map` is a single flat 1-based array: output global cell → the input
vertex it was built around. See `doc/tet_poly.md`.
"""
function tet_to_poly(m::Mesh; feature_angle::Real=1.0, point_data_to_cell_data::Bool=false,
                     record_vertex_ids::Bool=false)
    opts = _CTetToPolyOpts(Float64(feature_angle), Int32(point_data_to_cell_data),
                           Int32(record_vertex_ids), ntuple(_ -> Int64(0), 6))
    result = _check_ptr(ccall(_sym(:mio_tet_to_poly), Ptr{Cvoid},
                              (Ptr{Cvoid}, Ref{_CTetToPolyOpts}), _handle(m), Ref(opts)))
    try
        vm = _result_map(result, :mio_tet_to_poly_result_vertex_map)
        nnew, nbnd, nsplit = Ref{Int64}(0), Ref{Int64}(0), Ref{Int64}(0)
        _check(ccall(_sym(:mio_tet_to_poly_result_counts), Cint,
                     (Ptr{Cvoid}, Ptr{Int64}, Ptr{Int64}, Ptr{Int64}), result, nnew, nbnd, nsplit))
        out = Mesh(_check_ptr(ccall(_sym(:mio_tet_to_poly_result_take_mesh), Ptr{Cvoid},
                                    (Ptr{Cvoid},), result)))
        (mesh=out, vertex_map=vm, num_new_points=Int(nnew[]),
         num_boundary_faces=Int(nbnd[]), num_split_vertices=Int(nsplit[]))
    finally
        ccall(_sym(:mio_tet_to_poly_result_free), Cvoid, (Ptr{Cvoid},), result)
    end
end

"""
    poly_to_tet(mesh; record_parent_ids=false) -> (; mesh, point_map, cell_maps,
                                                   num_exact_tets, num_fanned_cells)

Decompose a polyhedral (or mixed) mesh into tetrahedra:
[`convert_cells`](@ref) with `"simplexify"`, except that a polyhedron that
already is a tetrahedron becomes that one tetra, with no new points. The maps
have `convert_cells`' shape. See `doc/tet_poly.md`.
"""
function poly_to_tet(m::Mesh; record_parent_ids::Bool=false)
    nexact, nfanned = Ref{Int64}(0), Ref{Int64}(0)
    result = _check_ptr(ccall(_sym(:mio_poly_to_tet), Ptr{Cvoid},
                              (Ptr{Cvoid}, Cint, Ptr{Int64}, Ptr{Int64}), _handle(m),
                              record_parent_ids ? Cint(1) : Cint(0), nexact, nfanned))
    try
        pm = _result_map(result, :mio_convert_cells_result_point_map)
        cm = _result_cell_maps(result, :mio_convert_cells_result_num_cell_maps,
                               :mio_convert_cells_result_cell_map)
        out = Mesh(_check_ptr(ccall(_sym(:mio_convert_cells_result_take_mesh), Ptr{Cvoid},
                                    (Ptr{Cvoid},), result)))
        (mesh=out, point_map=pm, cell_maps=cm, num_exact_tets=Int(nexact[]),
         num_fanned_cells=Int(nfanned[]))
    finally
        ccall(_sym(:mio_convert_cells_result_free), Cvoid, (Ptr{Cvoid},), result)
    end
end

"""
    blend_steps(a, b, w; blend_points=false) -> Mesh

`a` with every floating-point data array replaced by `(1 - w) a + w b`; the two
steps must share a topology. The kernel sequence resampling builds on. See
`doc/sequences.md`.
"""
function blend_steps(a::Mesh, b::Mesh, w::Real; blend_points::Bool=false)
    Mesh(_check_ptr(ccall(_sym(:mio_blend_steps), Ptr{Cvoid},
                          (Ptr{Cvoid}, Ptr{Cvoid}, Cdouble, Int32),
                          _handle(a), _handle(b), Float64(w), Int32(blend_points))))
end

const _REFINE_CLOSURES = Dict("" => Int32(0), "redgreen" => Int32(0), "red-green" => Int32(0),
                              "green" => Int32(0), "propagate" => Int32(1), "red" => Int32(1),
                              "balanced" => Int32(2), "2:1" => Int32(2))

const _REFINE_COMPARES = Dict("<" => Int32(0), "lt" => Int32(0), "<=" => Int32(1),
                              "le" => Int32(1), ">" => Int32(2), "gt" => Int32(2),
                              ">=" => Int32(3), "ge" => Int32(3), "==" => Int32(4),
                              "=" => Int32(4), "eq" => Int32(4), "!=" => Int32(5),
                              "ne" => Int32(5))

"""
    refine(mesh; levels=1, record_parent_ids=false, cells=nothing, region="",
           where_array="", where_op="<", where_value=0.0, closure="redgreen",
           record_levels=false, record_hierarchy=false)
        -> (; mesh, point_map, cell_maps)

Subdivide cells into same-type children (line → 2, triangle → 4, quad → 4,
tetra → 8, wedge → 8, hexahedron → 8), sharing the new edge/face/body nodes so
the result has no hanging nodes.

With no selector every cell is refined. Give at most **one** of `cells` (global,
block-major, **1-based** indices), `region` (a cell region selects its own cells;
a point region every cell with any node in it; a side region is an error) or
`where_array` + `where_op` + `where_value` (a threshold on a scalar `cell_data`
array), and only those cells are refined — the hanging nodes that leaves are
resolved by `closure`, so the output is still conforming. `"redgreen"` keeps the
extra refinement local; `"propagate"` is always available for every cell type but
converges to uniform refinement of the whole edge-connected component;
`"balanced"` keeps the hanging nodes and only enforces 2:1 balance, so the output
is **not conforming** but the cost is bounded by the selection (the constrained
nodes come back in `refine:hanging`).

Higher-order cells, pyramids and ragged blocks have no same-type subdivision
and fail by name. See `doc/refine.md`.

`record_hierarchy` attaches the `refine:cell_id`/`refine:parent_id` `cell_data`
arrays -- the persistent parent/child hierarchy a multigrid caller resolves
across the sequence of meshes it keeps ("a link between two meshes, not a tree
inside one"): an unsplit cell keeps its id and is its own parent; a split
cell's children each get a fresh id and carry the parent's id. An input
already carrying `refine:cell_id` is updated whatever this says. Also forces
`refine:entity` to be attached even when the closure leaves no hanging node,
since it already records the coarse corners each new fine node is the mean of
-- the multigrid prolongation weights, which `"redgreen"`/`"propagate"` would
otherwise never expose.
"""
function refine(m::Mesh; levels::Integer=1, record_parent_ids::Bool=false,
                cells=nothing, region::AbstractString="",
                where_array::AbstractString="", where_op::AbstractString="<",
                where_value::Real=0.0, closure::AbstractString="redgreen",
                record_levels::Bool=false, record_hierarchy::Bool=false)
    closure_code = get(_REFINE_CLOSURES, String(closure)) do
        throw(ArgumentError("refine: unknown closure '$closure' " *
                            "(expected 'redgreen'/'green', 'propagate' or 'balanced')"))
    end
    compare_code = get(_REFINE_COMPARES, String(where_op)) do
        throw(ArgumentError("refine: unknown comparison '$where_op' " *
                            "(expected '<', '<=', '>', '>=', '==' or '!=')"))
    end
    # 1-based on this side, 0-based across the ABI -- the binding's standing rule
    # for every index array.
    ids = cells === nothing ? Int64[] : Vector{Int64}(vec(collect(cells))) .- Int64(1)
    region_c = Vector{UInt8}(codeunits(String(region)))
    push!(region_c, 0x00)
    array_c = Vector{UInt8}(codeunits(String(where_array)))
    push!(array_c, 0x00)
    result = GC.@preserve ids region_c array_c begin
        opts = _CRefineOpts(isempty(ids) ? Ptr{Int64}(C_NULL) : pointer(ids),
                            Int64(length(ids)),
                            Cstring(pointer(region_c)), Cstring(pointer(array_c)),
                            Cdouble(where_value), Int32(levels),
                            record_parent_ids ? Int32(1) : Int32(0),
                            record_levels ? Int32(1) : Int32(0),
                            closure_code, compare_code, Int32(0),
                            record_hierarchy ? Int64(1) : Int64(0),
                            (Int64(0), Int64(0), Int64(0), Int64(0), Int64(0)))
        _check_ptr(ccall(_sym(:mio_refine_ex), Ptr{Cvoid},
                         (Ptr{Cvoid}, Ref{_CRefineOpts}), _handle(m), Ref(opts)))
    end
    try
        pm = _result_map(result, :mio_refine_result_point_map)
        cm = _result_cell_maps(result, :mio_refine_result_num_cell_maps,
                               :mio_refine_result_cell_map)
        out = Mesh(_check_ptr(ccall(_sym(:mio_refine_result_take_mesh), Ptr{Cvoid},
                                    (Ptr{Cvoid},), result)))
        (mesh=out, point_map=pm, cell_maps=cm)
    finally
        ccall(_sym(:mio_refine_result_free), Cvoid, (Ptr{Cvoid},), result)
    end
end

"""
    tessellate(m; levels=2, curved=true, fields=true, record_stencil=false)
        -> (; mesh, num_source_points, num_source_cells, num_curved_source_cells,
             num_pass_through_source_cells, num_points, num_cells,
             num_boundary_facets, num_interior_facets, num_facets_with_bad_count)

Tessellate the curved (quadratic) cells isoparametrically: every `triangle6`,
`quad8`, `quad9`, `tetra10` or `hexahedron27` cell is sampled on a reference
lattice and mapped through its own shape functions (curved surface cells become
triangles, curved volume cells tetras); every other cell passes through. The
output carries point data `tessellate:source_point` (the 0-based source point,
or -1) and cell data `tessellate:source_cell` / `tessellate:sub_index`;
`record_stencil=true` adds the `(27, num_points)` `tessellate:stencil` and
`tessellate:weights`. `num_facets_with_bad_count == 0` for a conforming result.
See `doc/tessellation.md`.
"""
function tessellate(m::Mesh; levels::Integer=2, curved::Bool=true, fields::Bool=true,
                    record_stencil::Bool=false)
    report = Ref{_CTessellateReport}()
    opts = _CTessellateOpts(Int64(levels), Int32(curved), Int32(fields),
                            Int32(record_stencil), Int32(0), ntuple(_ -> Int64(0), 6))
    ptr = ccall(_sym(:mio_tessellate), Ptr{Cvoid},
                (Ptr{Cvoid}, Ref{_CTessellateOpts}, Ptr{_CTessellateReport}),
                _handle(m), Ref(opts), report)
    r = _check_ptr(ptr)
    rep = report[]
    (mesh=Mesh(r), num_source_points=Int(rep.num_source_points),
     num_source_cells=Int(rep.num_source_cells),
     num_curved_source_cells=Int(rep.num_curved_source_cells),
     num_pass_through_source_cells=Int(rep.num_pass_through_source_cells),
     num_points=Int(rep.num_points), num_cells=Int(rep.num_cells),
     num_boundary_facets=Int(rep.num_boundary_facets),
     num_interior_facets=Int(rep.num_interior_facets),
     num_facets_with_bad_count=Int(rep.num_facets_with_bad_count))
end

"""
    decimate(mesh; ratio=-1.0, target_faces=-1, max_error=-1.0, placement="optimal",
             preserve_boundary=true, preserve_features=true, feature_angle=30.0,
             max_deviation=-1.0, coverage=0.999, max_faces=-1)
        -> (; mesh, point_map, cell_maps, faces_removed, points_removed,
              collapses_rejected, max_error_applied, deviation, deviation_met, rounds)

Reduce a **surface** mesh's face count by greedy quadric-error-metric edge
collapse — the inverse of [`refine`](@ref). Exactly one of `ratio` (fraction of
faces to KEEP), `target_faces`, `max_error` and `max_deviation` (the largest
two-sided distance from the input on a `coverage` share of samples, roadmap
7.8.5) must be set; a negative value means "unset". `max_faces` caps any of
them (alone, it is the target).

A 3-D volume mesh fails by name (extract the surface first), as do
higher-order, ragged and line/vertex blocks. `frozen` is an optional vector
of **1-based** point ids to pin outright; shifted to 0-based across the C
ABI, as every index array in this binding. An out-of-range id raises by
name. See `doc/decimate.md`.
"""
function decimate(m::Mesh; ratio::Real=-1.0, target_faces::Integer=-1,
                  max_error::Real=-1.0, placement::AbstractString="optimal",
                  preserve_boundary::Bool=true, preserve_features::Bool=true,
                  feature_angle::Real=30.0, frozen=nothing, max_deviation::Real=-1.0,
                  coverage::Real=0.999, max_faces::Integer=-1)
    placement_c = Vector{UInt8}(codeunits(String(placement)))
    push!(placement_c, 0x00)
    ids = frozen === nothing ? Int64[] : Vector{Int64}(vec(collect(frozen))) .- Int64(1)
    result = GC.@preserve ids placement_c begin
        opts = _CDecimateOpts(Cdouble(ratio), Int64(target_faces), Cdouble(max_error),
                              Cstring(pointer(placement_c)),
                              preserve_boundary ? Int32(1) : Int32(0),
                              preserve_features ? Int32(1) : Int32(0),
                              Cdouble(feature_angle),
                              isempty(ids) ? Ptr{Int64}(C_NULL) : pointer(ids),
                              Int64(length(ids)),
                              max_deviation >= 0 ? Cdouble(max_deviation) : Cdouble(0),
                              Cdouble(coverage), max_faces > 0 ? Int64(max_faces) : Int64(0),
                              (Int64(0),))
        _check_ptr(ccall(_sym(:mio_decimate_ex), Ptr{Cvoid},
                         (Ptr{Cvoid}, Ref{_CDecimateOpts}), _handle(m), Ref(opts)))
    end
    try
        pm = _result_map(result, :mio_decimate_result_point_map)
        cm = _result_cell_maps(result, :mio_decimate_result_num_cell_maps,
                               :mio_decimate_result_cell_map)
        fr = ccall(_sym(:mio_decimate_result_faces_removed), Int64, (Ptr{Cvoid},), result)
        pr = ccall(_sym(:mio_decimate_result_points_removed), Int64, (Ptr{Cvoid},), result)
        cr = ccall(_sym(:mio_decimate_result_collapses_rejected), Int64, (Ptr{Cvoid},), result)
        me = ccall(_sym(:mio_decimate_result_max_error_applied), Cdouble, (Ptr{Cvoid},), result)
        met = Ref{Int32}(0)
        rounds = Ref{Int64}(0)
        dev = ccall(_sym(:mio_decimate_result_deviation), Cdouble,
                    (Ptr{Cvoid}, Ref{Int32}, Ref{Int64}), result, met, rounds)
        out = Mesh(_check_ptr(ccall(_sym(:mio_decimate_result_take_mesh), Ptr{Cvoid},
                                    (Ptr{Cvoid},), result)))
        (mesh=out, point_map=pm, cell_maps=cm, faces_removed=Int(fr),
         points_removed=Int(pr), collapses_rejected=Int(cr),
         max_error_applied=Float64(me), deviation=Float64(dev), deviation_met=met[] != 0,
         rounds=Int(rounds[]))
    finally
        ccall(_sym(:mio_decimate_result_free), Cvoid, (Ptr{Cvoid},), result)
    end
end

"""
    decimate_volume(mesh; ratio=-1.0, target_cells=-1, max_error=-1.0, placement="optimal",
                    preserve_boundary=false, preserve_features=true, feature_angle=30.0,
                    frozen=nothing)
        -> (; mesh, point_map, cell_maps, tets_removed, points_removed,
              collapses_rejected, max_error_applied)

Reduce a **tet** mesh's cell count by greedy quadric-error-metric tet-edge
collapse — the volume-mesh sibling of [`decimate`](@ref). Exactly one of
`ratio` (fraction of tets to KEEP), `target_cells` and `max_error` must be set;
a negative value means "unset".

The mesh must be tet-only: hexahedron/wedge/pyramid/polyhedron, ragged,
higher-order or mixed lower-dimensional blocks fail by name (simplexify first
with [`convert_cells`](@ref)). Unlike `decimate`, boundary vertices participate
by default (`preserve_boundary=false`). `frozen` is an optional vector of
**1-based** point ids to pin outright; shifted to 0-based across the C ABI, as
every index array in this binding. An out-of-range id raises by name. The maps
are 1-based with 0 for "absent". See `doc/decimate_volume.md`.
"""
function decimate_volume(m::Mesh; ratio::Real=-1.0, target_cells::Integer=-1,
                         max_error::Real=-1.0, placement::AbstractString="optimal",
                         preserve_boundary::Bool=false, preserve_features::Bool=true,
                         feature_angle::Real=30.0, frozen=nothing)
    placement_c = Vector{UInt8}(codeunits(String(placement)))
    push!(placement_c, 0x00)
    ids = frozen === nothing ? Int64[] : Vector{Int64}(vec(collect(frozen))) .- Int64(1)
    result = GC.@preserve ids placement_c begin
        opts = _CDecimateVolumeOpts(Cdouble(ratio), Int64(target_cells), Cdouble(max_error),
                                    Cstring(pointer(placement_c)),
                                    preserve_boundary ? Int32(1) : Int32(0),
                                    preserve_features ? Int32(1) : Int32(0),
                                    Cdouble(feature_angle),
                                    isempty(ids) ? Ptr{Int64}(C_NULL) : pointer(ids),
                                    Int64(length(ids)), (Int64(0), Int64(0), Int64(0), Int64(0)))
        _check_ptr(ccall(_sym(:mio_decimate_volume_ex), Ptr{Cvoid},
                         (Ptr{Cvoid}, Ref{_CDecimateVolumeOpts}), _handle(m), Ref(opts)))
    end
    try
        pm = _result_map(result, :mio_decimate_volume_result_point_map)
        cm = _result_cell_maps(result, :mio_decimate_volume_result_num_cell_maps,
                               :mio_decimate_volume_result_cell_map)
        tr = ccall(_sym(:mio_decimate_volume_result_tets_removed), Int64, (Ptr{Cvoid},), result)
        pr = ccall(_sym(:mio_decimate_volume_result_points_removed), Int64, (Ptr{Cvoid},),
                   result)
        cr = ccall(_sym(:mio_decimate_volume_result_collapses_rejected), Int64, (Ptr{Cvoid},),
                   result)
        me = ccall(_sym(:mio_decimate_volume_result_max_error_applied), Cdouble, (Ptr{Cvoid},),
                   result)
        out = Mesh(_check_ptr(ccall(_sym(:mio_decimate_volume_result_take_mesh), Ptr{Cvoid},
                                    (Ptr{Cvoid},), result)))
        (mesh=out, point_map=pm, cell_maps=cm, tets_removed=Int(tr),
         points_removed=Int(pr), collapses_rejected=Int(cr),
         max_error_applied=Float64(me))
    finally
        ccall(_sym(:mio_decimate_volume_result_free), Cvoid, (Ptr{Cvoid},), result)
    end
end

"""
    partition(mesh, nparts; method="auto", imbalance=0.03, mode="eco", seed=0,
              record_ids=false, ghost_layers=0, weights="")
        -> Vector{@NamedTuple{part_id::Int, mesh::Mesh, point_map::Vector{Int},
                              cell_maps::Vector{Vector{Int}}}}

Decompose a mesh into exactly `nparts` balanced pieces (the count-driven
complement to [`split`](@ref)). `method` is `"sfc"` (Hilbert curve, always
available), `"kahip"` (only in a build configured with KaHIP; requesting it
elsewhere fails by name) or `"auto"`.

Every piece keeps the input's cell-block structure 1:1, so the cell maps index
by input block and the pieces recombine into the input. `ghost_layers > 0`
instead grows each piece by that many shared-node BFS layers of other parts'
cells (a halo), tagged `partition:ghost` (0 = owned, L = reached at layer L),
in which case the pieces overlap. See `doc/partition.md`.
"""
function partition(m::Mesh, nparts::Integer; method::AbstractString="auto",
                   imbalance::Real=0.03, mode::AbstractString="eco", seed::Integer=0,
                   record_ids::Bool=false, ghost_layers::Integer=0,
                   weights::AbstractString="")
    result = _check_ptr(ccall(_sym(:mio_partition), Ptr{Cvoid},
                              (Ptr{Cvoid}, Cint, Cstring, Cdouble, Cstring, Cint, Cint,
                               Cint, Cstring),
                              _handle(m), Cint(nparts), method, Float64(imbalance), mode,
                              Cint(seed), record_ids ? Cint(1) : Cint(0),
                              Cint(ghost_layers), weights))
    try
        n = _check_count(ccall(_sym(:mio_partition_result_num_pieces), Int64,
                               (Ptr{Cvoid},), result), "mio_partition_result_num_pieces")
        nmaps = _check_count(ccall(_sym(:mio_partition_result_num_cell_maps), Int64,
                                   (Ptr{Cvoid},), result),
                             "mio_partition_result_num_cell_maps")
        out = @NamedTuple{part_id::Int, mesh::Mesh, point_map::Vector{Int},
                          cell_maps::Vector{Vector{Int}}}[]
        for i in 0:(n-1)
            pid = ccall(_sym(:mio_partition_result_part_id), Cint, (Ptr{Cvoid}, Int64),
                        result, i)
            data = Ref{Ptr{Cvoid}}(C_NULL); dt = Ref{Cint}(0); len = Ref{Int64}(0)
            _check(ccall(_sym(:mio_partition_result_point_map), Cint,
                         (Ptr{Cvoid}, Int64, Ptr{Ptr{Cvoid}}, Ptr{Cint}, Ptr{Int64}),
                         result, i, data, dt, len))
            pm = _shift_map(data[], len[])
            cms = Vector{Int}[]
            for b in 0:(nmaps-1)
                d2 = Ref{Ptr{Cvoid}}(C_NULL); t2 = Ref{Cint}(0); l2 = Ref{Int64}(0)
                _check(ccall(_sym(:mio_partition_result_cell_map), Cint,
                             (Ptr{Cvoid}, Int64, Int64, Ptr{Ptr{Cvoid}}, Ptr{Cint},
                              Ptr{Int64}), result, i, b, d2, t2, l2))
                push!(cms, _shift_map(d2[], l2[]))
            end
            piece = Mesh(_check_ptr(ccall(_sym(:mio_partition_result_take_mesh), Ptr{Cvoid},
                                          (Ptr{Cvoid}, Int64), result, i)))
            push!(out, (part_id=Int(pid), mesh=piece, point_map=pm, cell_maps=cms))
        end
        out
    finally
        ccall(_sym(:mio_partition_result_free), Cvoid, (Ptr{Cvoid},), result)
    end
end

"""
    partition_labels(mesh, nparts; method="auto", imbalance=0.03, mode="eco",
                     seed=0, weights="") -> Vector{Int}

The per-cell part assignment only, without building the pieces: one entry per
cell in **block-major global cell order** (block 1's cells first, then block
2's, …), with values in `0:nparts-1`.

These are part **ids**, not indices, so they are deliberately **not** shifted
to 1-based — the same rule the Fortran module follows.
"""
function partition_labels(m::Mesh, nparts::Integer; method::AbstractString="auto",
                          imbalance::Real=0.03, mode::AbstractString="eco",
                          seed::Integer=0, weights::AbstractString="")
    n = num_cells(m)
    labels = Vector{Int64}(undef, n)
    _check(GC.@preserve labels ccall(_sym(:mio_partition_labels), Cint,
                                     (Ptr{Cvoid}, Cint, Cstring, Cdouble, Cstring, Cint,
                                      Cstring, Ptr{Int64}, Int64),
                                     _handle(m), Cint(nparts), method, Float64(imbalance),
                                     mode, Cint(seed), weights, pointer(labels), Int64(n)))
    Int[Int(x) for x in labels]
end

# --- data operations ---------------------------------------------------------
#
# These act on point_data / cell_data / field_data and never touch geometry:
# points, connectivity and block order come through bit-identical.

const _LOCATIONS = Dict(:point => MIO_DATA_POINT, :cell => MIO_DATA_CELL,
                        :field => MIO_DATA_FIELD)

function _location(loc::Symbol)
    haskey(_LOCATIONS, loc) || throw(ArgumentError(
        "location must be :point, :cell or :field, got :$loc"))
    _LOCATIONS[loc]
end

"""
    data_drop(mesh, location, names; ignore_missing=false) -> Mesh

Drop the named data arrays at one location (`:point`, `:cell` or `:field`).

The C ABI exposes the three primitives (`drop` / `keep` / `rename`) but not the
combined `data_manage`; they compose to the same effect. See
`doc/data_manage.md`.
"""
function data_drop(m::Mesh, location::Symbol, names; ignore_missing::Bool=false)
    h = _handle(m); loc = _location(location)
    ptr = _with_names(names) do p, c
        ccall(_sym(:mio_data_drop), Ptr{Cvoid},
              (Ptr{Cvoid}, Cint, Ptr{Cstring}, Int64, Cint), h, loc, p, c,
              ignore_missing ? Cint(1) : Cint(0))
    end
    Mesh(_check_ptr(ptr))
end

"""
    sets_to_data(mesh, location; data_name=nothing, join_char="-", order=String[]) -> Mesh

Convert point/cell region-backed sets to scalar integer labels on a new mesh.
Labels are zero-based, later overlapping sets win, and uncovered rows are -1.
`order` must name every set once, or be empty to use native region-name order.
"""
function sets_to_data(m::Mesh, location::Symbol; data_name=nothing,
                      join_char::AbstractString="-", order=String[])
    h = _handle(m); loc = _location(location)
    ptr = _with_names(order) do p, c
        ccall(_sym(:mio_sets_to_data), Ptr{Cvoid},
              (Ptr{Cvoid}, Cint, Cstring, Cstring, Ptr{Cstring}, Int64),
              h, loc, isnothing(data_name) ? C_NULL : data_name, join_char, p, c)
    end
    Mesh(_check_ptr(ptr))
end

"""
    data_to_sets(mesh, location, key) -> Mesh

Convert a scalar integer point/cell field into region-backed sets, removing
the source field on the returned mesh. Geometry and unrelated metadata survive.
"""
data_to_sets(m::Mesh, location::Symbol, key::AbstractString) =
    Mesh(_check_ptr(ccall(_sym(:mio_data_to_sets), Ptr{Cvoid},
                          (Ptr{Cvoid}, Cint, Cstring), _handle(m), _location(location), key)))

"""
    data_keep(mesh, location, names; ignore_missing=false) -> Mesh

Keep only the named arrays at one location, dropping the rest **there**; the
other two locations are untouched. An empty `names` drops everything at
`location`.
"""
function data_keep(m::Mesh, location::Symbol, names; ignore_missing::Bool=false)
    h = _handle(m); loc = _location(location)
    ptr = _with_names(names) do p, c
        ccall(_sym(:mio_data_keep), Ptr{Cvoid},
              (Ptr{Cvoid}, Cint, Ptr{Cstring}, Int64, Cint), h, loc, p, c,
              ignore_missing ? Cint(1) : Cint(0))
    end
    Mesh(_check_ptr(ptr))
end

"""
    data_rename(mesh, location, from, to) -> Mesh

Rename one data array, preserving its values, dtype and shape.
"""
data_rename(m::Mesh, location::Symbol, from::AbstractString, to::AbstractString) =
    Mesh(_check_ptr(ccall(_sym(:mio_data_rename), Ptr{Cvoid},
                          (Ptr{Cvoid}, Cint, Cstring, Cstring),
                          _handle(m), _location(location), from, to)))

"""
    data_point_to_cell(mesh, names=String[]; suffix="") -> Mesh

Average point data onto the cells: each cell's value is the mean over its own
nodes. The output is always `Float64`. An empty `names` converts every
point-data array. See `doc/data_average.md`.
"""
function data_point_to_cell(m::Mesh, names=String[]; suffix::AbstractString="")
    h = _handle(m)
    ptr = _with_names(names) do p, c
        ccall(_sym(:mio_data_point_to_cell), Ptr{Cvoid},
              (Ptr{Cvoid}, Ptr{Cstring}, Int64, Cstring), h, p, c, suffix)
    end
    Mesh(_check_ptr(ptr))
end

"""
    data_cell_to_point(mesh, names=String[]; weight=:uniform, suffix="") -> Mesh

Average cell data onto the points: each point's value is the (optionally
measure-weighted) mean over the incident cells. `weight` is `:uniform` or
`:measure`. A point touched by no cell gets `NaN`; the output is always
`Float64`.
"""
function data_cell_to_point(m::Mesh, names=String[]; weight::Symbol=:uniform,
                            suffix::AbstractString="")
    w = weight === :uniform ? Cint(0) : weight === :measure ? Cint(1) :
        throw(ArgumentError("weight must be :uniform or :measure, got :$weight"))
    h = _handle(m)
    ptr = _with_names(names) do p, c
        ccall(_sym(:mio_data_cell_to_point), Ptr{Cvoid},
              (Ptr{Cvoid}, Ptr{Cstring}, Int64, Cint, Cstring), h, p, c, w, suffix)
    end
    Mesh(_check_ptr(ptr))
end

"""
    data_calc(mesh, expression, output_name; location=:point, overwrite=false) -> Mesh

Evaluate an elementwise expression over the arrays at `location` and store the
result there as `output_name`.

The grammar accepts `+ - * /`, unary minus, parentheses, numeric literals,
array names, and the functions `abs`, `sqrt`, `min`, `max` and `norm`. Nothing
else is evaluated — there is **no arbitrary-code path**. See `doc/data_calc.md`.
"""
data_calc(m::Mesh, expression::AbstractString, output_name::AbstractString;
          location::Symbol=:point, overwrite::Bool=false) =
    Mesh(_check_ptr(ccall(_sym(:mio_data_calc), Ptr{Cvoid},
                          (Ptr{Cvoid}, Cstring, Cint, Cstring, Cint),
                          _handle(m), expression, _location(location), output_name,
                          overwrite ? Cint(1) : Cint(0))))

"""
    data_condition(mesh, location, names=String[]; mode=:clamp, lo=0.0, hi=1.0,
                   scope=:component, nan_policy=:ignore, nan_replacement=0.0,
                   suffix="") -> Mesh

Condition the values of the selected arrays. `mode` is `:clamp`, `:normalize`
or `:standardize`; `scope` is `:component` (each trailing component
independently) or `:magnitude` (rescale whole rows, preserving direction);
`nan_policy` is `:ignore`, `:replace` or `:fail`.

For cell data the statistics are computed **jointly across all blocks**. An
empty `names` conditions every array at `location`. See
`doc/data_condition.md`.
"""
function data_condition(m::Mesh, location::Symbol, names=String[]; mode::Symbol=:clamp,
                        lo::Real=0.0, hi::Real=1.0, scope::Symbol=:component,
                        nan_policy::Symbol=:ignore, nan_replacement::Real=0.0,
                        suffix::AbstractString="")
    md = mode === :clamp ? Cint(0) : mode === :normalize ? Cint(1) :
         mode === :standardize ? Cint(2) :
         throw(ArgumentError("mode must be :clamp, :normalize or :standardize"))
    sc = scope === :component ? Cint(0) : scope === :magnitude ? Cint(1) :
         throw(ArgumentError("scope must be :component or :magnitude"))
    np = nan_policy === :ignore ? Cint(0) : nan_policy === :replace ? Cint(1) :
         nan_policy === :fail ? Cint(2) :
         throw(ArgumentError("nan_policy must be :ignore, :replace or :fail"))
    h = _handle(m); loc = _location(location)
    ptr = _with_names(names) do p, c
        ccall(_sym(:mio_data_condition), Ptr{Cvoid},
              (Ptr{Cvoid}, Cint, Ptr{Cstring}, Int64, Cint, Cdouble, Cdouble, Cint, Cint,
               Cdouble, Cstring),
              h, loc, p, c, md, Float64(lo), Float64(hi), sc, np,
              Float64(nan_replacement), suffix)
    end
    Mesh(_check_ptr(ptr))
end

"""
    tensor_invariants(mesh, location, names=String[]; outputs=Symbol[], prefix="",
                      suffix="", overwrite=true) -> Mesh

von Mises, principal, hydrostatic and deviatoric fields of a symmetric
(6-component, `xx yy zz xy yz zx`) or general 3x3 (9-component, row-major)
tensor array. `mises`/`principal` use the symmetric part of a 9-component
input. An empty `names` processes every 6- or 9-component array at `location`
(`:field` is rejected). `outputs` is any of `:mises`, `:principal`,
`:hydrostatic`, `:deviatoric`; an empty vector means all four. Results are
stored as `prefix * name * "_" * output * suffix`. See
`doc/tensor_invariants.md`.
"""
function tensor_invariants(m::Mesh, location::Symbol, names=String[]; outputs=Symbol[],
                           prefix::AbstractString="", suffix::AbstractString="",
                           overwrite::Bool=true)
    flags = Dict(:mises => 1, :principal => 2, :hydrostatic => 4, :deviatoric => 8)
    mask = isempty(outputs) ? 15 :
           reduce(|, (get(flags, o) do
                          throw(ArgumentError("outputs entries must be :mises, :principal, " *
                                              ":hydrostatic or :deviatoric, got :$o"))
                      end for o in outputs))
    h = _handle(m); loc = _location(location)
    ptr = _with_names(names) do p, c
        ccall(_sym(:mio_tensor_invariants), Ptr{Cvoid},
              (Ptr{Cvoid}, Cint, Ptr{Cstring}, Int64, Cint, Cstring, Cstring, Cint),
              h, loc, p, c, Cint(mask), prefix, suffix, overwrite ? Cint(1) : Cint(0))
    end
    Mesh(_check_ptr(ptr))
end

"""
    data_info(mesh) -> Vector{NamedTuple}

Summarize every data array the mesh carries: location, name, dtype, shape,
counts, whole-array and per-component min/max/mean, and the NaN/inf/finite
counts. Read-only, and it **never throws on non-finite values** — it counts
them. See `doc/data_info.md`.
"""
function data_info(m::Mesh)
    handle = _check_ptr(ccall(_sym(:mio_data_info_create), Ptr{Cvoid}, (Ptr{Cvoid},),
                              _handle(m)))
    try
        n = _check_count(ccall(_sym(:mio_data_info_count), Int64, (Ptr{Cvoid},), handle),
                         "mio_data_info_count")
        out = []
        for i in 0:(n-1)
            name = _getstring() do buf, len
                ccall(_sym(:mio_data_info_name), Int64,
                      (Ptr{Cvoid}, Int64, Ptr{UInt8}, Int64), handle, i, buf, len)
            end
            entry = Ref{_CDataArrayInfo}()
            _check(ccall(_sym(:mio_data_info_entry), Cint,
                         (Ptr{Cvoid}, Int64, Ptr{_CDataArrayInfo}), handle, i, entry))
            e = entry[]
            comps = NamedTuple{(:min, :max, :mean),Tuple{Float64,Float64,Float64}}[]
            for c in 0:(e.num_components-1)
                lo = Ref{Cdouble}(0.0); hi = Ref{Cdouble}(0.0); mu = Ref{Cdouble}(0.0)
                _check(ccall(_sym(:mio_data_info_component), Cint,
                             (Ptr{Cvoid}, Int64, Int64, Ptr{Cdouble}, Ptr{Cdouble},
                              Ptr{Cdouble}), handle, i, c, lo, hi, mu))
                push!(comps, (min=Float64(lo[]), max=Float64(hi[]), mean=Float64(mu[])))
            end
            push!(out, (name=name, location=_LOCATION_SYMS[e.location+1],
                        dtype=_jltype(e.dtype), num_blocks=Int(e.num_blocks),
                        num_entries=Int(e.num_entries),
                        num_components=Int(e.num_components),
                        num_values=Int(e.num_values), min=e.min, max=e.max, mean=e.mean,
                        num_nan=Int(e.num_nan), num_inf=Int(e.num_inf),
                        num_finite=Int(e.num_finite),
                        inconsistent_blocks=e.inconsistent_blocks != 0,
                        components=comps))
        end
        out
    finally
        ccall(_sym(:mio_data_info_free), Cvoid, (Ptr{Cvoid},), handle)
    end
end

function _field_integral_components(handle, i::Int64, region, nc::Int64)
    comps = NamedTuple{(:total, :mean, :domain_measure, :num_nan),
                        Tuple{Float64,Float64,Float64,Int64}}[]
    for c in 0:(nc - 1)
        total = Ref{Cdouble}(0.0); mean = Ref{Cdouble}(0.0)
        dm = Ref{Cdouble}(0.0); nnan = Ref{Int64}(0)
        if region === nothing
            _check(ccall(_sym(:mio_data_integrate_component), Cint,
                         (Ptr{Cvoid}, Int64, Int64, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble},
                          Ptr{Int64}), handle, i, c, total, mean, dm, nnan))
        else
            _check(ccall(_sym(:mio_data_integrate_region_component), Cint,
                         (Ptr{Cvoid}, Int64, Int64, Int64, Ptr{Cdouble}, Ptr{Cdouble},
                          Ptr{Cdouble}, Ptr{Int64}), handle, i, region, c, total, mean, dm, nnan))
        end
        push!(comps, (total=total[], mean=mean[], domain_measure=dm[], num_nan=Int(nnan[])))
    end
    comps
end

"""
    data_integrate(mesh, names=String[]) -> Vector{NamedTuple}

Cell-measure-weighted total/mean of one or more `cell_data` arrays --
`gradient`'s integration counterpart (`gradient` differentiates a field, this
integrates one). Every sum is weighted by the cell's own length/area/volume,
excluding both a cell with unmeasurable geometry and a component with a
non-finite value from that component's numerator *and* denominator. `names`
empty means every `cell_data` array. Reported for the whole mesh (`domain`)
and independently for every named Cell region (`regions`) -- a cell in two
regions contributes fully to both, one in none contributes to neither. A
`point_data`-only name throws naming `data_point_to_cell` as the fix. See
`doc/field_integration.md`.
"""
function data_integrate(m::Mesh, names=String[])
    h = _handle(m)
    handle = _with_names(names) do p, c
        _check_ptr(ccall(_sym(:mio_data_integrate_create), Ptr{Cvoid},
                          (Ptr{Cvoid}, Ptr{Cstring}, Int64), h, p, c))
    end
    try
        n = _check_count(ccall(_sym(:mio_data_integrate_count), Int64, (Ptr{Cvoid},), handle),
                         "mio_data_integrate_count")
        out = []
        for i in 0:(n - 1)
            name = _getstring() do buf, len
                ccall(_sym(:mio_data_integrate_name), Int64,
                      (Ptr{Cvoid}, Int64, Ptr{UInt8}, Int64), handle, i, buf, len)
            end
            entry = Ref{_CFieldIntegralInfo}()
            _check(ccall(_sym(:mio_data_integrate_entry), Cint,
                         (Ptr{Cvoid}, Int64, Ptr{_CFieldIntegralInfo}), handle, i, entry))
            e = entry[]
            domain = (num_cells=Int(e.num_cells), num_skipped=Int(e.num_skipped),
                      components=_field_integral_components(handle, i, nothing, e.num_components))

            nregions = _check_count(ccall(_sym(:mio_data_integrate_region_count), Int64,
                                          (Ptr{Cvoid}, Int64), handle, i),
                                    "mio_data_integrate_region_count")
            regions = []
            for r in 0:(nregions - 1)
                rname = _getstring() do buf, len
                    ccall(_sym(:mio_data_integrate_region_name), Int64,
                          (Ptr{Cvoid}, Int64, Int64, Ptr{UInt8}, Int64), handle, i, r, buf, len)
                end
                rentry = Ref{_CFieldIntegralInfo}()
                _check(ccall(_sym(:mio_data_integrate_region_entry), Cint,
                             (Ptr{Cvoid}, Int64, Int64, Ptr{_CFieldIntegralInfo}), handle, i, r,
                             rentry))
                re = rentry[]
                push!(regions,
                      (name=rname, num_cells=Int(re.num_cells), num_skipped=Int(re.num_skipped),
                       components=_field_integral_components(handle, i, r, re.num_components)))
            end
            push!(out, (name=name, num_components=Int(e.num_components), domain=domain,
                        regions=regions))
        end
        out
    finally
        ccall(_sym(:mio_data_integrate_free), Cvoid, (Ptr{Cvoid},), handle)
    end
end

const _LOCATION_SYMS = (:point, :cell, :field)

# ---------------------------------------------------------------------------
# Settings pipeline (v9.11.0): run a whole settings.json (read -> operation
# chain -> write; PascalCase vocabulary, see doc/pipeline.md). The flat ABI
# carries JSON text only.
# ---------------------------------------------------------------------------

"""
    run_pipeline_file(settings_path)

Run the settings pipeline stored in the file at `settings_path` (a
`settings.json`; PascalCase ops/keys, see `doc/pipeline.md`). Needs a build
with the JSON parser (`-DMESHIOPLUSPLUS_WITH_JSON=ON`); otherwise the raised
error names the flag — [`pipeline_has_json`](@ref) reports which build this is.
"""
function run_pipeline_file(settings_path::AbstractString)
    _check(ccall(_sym(:mio_pipeline_run_file), Cint, (Cstring,), settings_path))
    nothing
end

"""
    run_pipeline_json(json_text)

Run a settings pipeline given as JSON text (see [`run_pipeline_file`](@ref)).
"""
function run_pipeline_json(json_text::AbstractString)
    _check(ccall(_sym(:mio_pipeline_run_json), Cint, (Cstring,), json_text))
    nothing
end

"""
    pipeline_has_json() -> Bool

Whether the loaded library carries the JSON pipeline parser.
"""
pipeline_has_json() = ccall(_sym(:mio_pipeline_has_json), Cint, ()) != 0

function _pipeline_report(run::Symbol, text::AbstractString)
    h = _check_ptr(ccall(_sym(run), Ptr{Cvoid}, (Cstring,), text))
    try
        _getstring((buf, len) -> ccall(_sym(:mio_pipeline_report_json), Int64,
                   (Ptr{Cvoid}, Ptr{UInt8}, Int64), h, buf, len))
    finally
        ccall(_sym(:mio_pipeline_report_free), Cvoid, (Ptr{Cvoid},), h)
    end
end

"""Run once and return the shared steps/warnings report as JSON text."""
run_pipeline_file_report(path::AbstractString) = _pipeline_report(:mio_pipeline_run_file_report, path)
"""Run JSON settings once and return steps/warnings as JSON text."""
run_pipeline_json_report(text::AbstractString) = _pipeline_report(:mio_pipeline_run_json_report, text)
"""Run a sequence settings file once and return steps/warnings as JSON text."""
run_sequence_file_report(path::AbstractString) = _pipeline_report(:mio_sequence_pipeline_run_file_report, path)
"""Run sequence JSON settings once and return steps/warnings as JSON text."""
run_sequence_json_report(text::AbstractString) = _pipeline_report(:mio_sequence_pipeline_run_json_report, text)

# --- regular grids and signed distance ---------------------------------------

const _SDF_SIGNS = Dict("unsigned" => Int32(0), "pseudonormal" => Int32(1),
                        "winding-number" => Int32(2), "winding_number" => Int32(2))
const _SDF_WEIGHTS = Dict("angle" => Int32(0), "area" => Int32(1))
const _SDF_LOCATIONS = Dict("corner" => Int32(0), "point" => Int32(0),
                            "center" => Int32(1), "centre" => Int32(1), "cell" => Int32(1))
const _SDF_CHECKS = Dict("off" => Int32(0), "warn" => Int32(1), "error" => Int32(2))
const _VOXEL_FILLS = Dict("all" => Int32(0), "surface" => Int32(1), "inside" => Int32(2))

function _code(table, value, what)
    key = replace(String(value), '_' => '-')
    haskey(table, key) && return table[key]
    haskey(table, String(value)) && return table[String(value)]
    throw(ArgumentError("meshio++: unknown $what '$(value)'"))
end

"""
    grid(dims; origin=(0,0,0), spacing=(1,1,1), max_cells=20_000_000) -> Mesh

Build a regular hexahedron lattice from nothing -- the only constructor here
that takes no input mesh. Points run x fastest, then y, then z.

A zero cell count on any axis gives an empty mesh, which is a legal request
rather than an error.
"""
function grid(dims; origin=(0.0, 0.0, 0.0), spacing=(1.0, 1.0, 1.0),
              max_cells::Integer=20_000_000)
    d = Int64[Int64(v) for v in dims]
    length(d) == 3 || throw(ArgumentError("dims must have three entries"))
    o = Cdouble[Cdouble(v) for v in origin]
    sp = Cdouble[Cdouble(v) for v in spacing]
    ptr = GC.@preserve d o sp ccall(_sym(:mio_grid), Ptr{Cvoid},
                                    (Ptr{Int64}, Ptr{Cdouble}, Ptr{Cdouble}, Int64),
                                    pointer(d), pointer(o), pointer(sp), Int64(max_cells))
    Mesh(_check_ptr(ptr))
end

# --- parametric generators ----------------------------------------------------

const _GEN_CELLS = Dict("" => Int32(0), "quad" => Int32(1), "triangle" => Int32(2),
                        "hexahedron" => Int32(3), "tetra" => Int32(4), "wedge" => Int32(5))
const _GEN_SPLITS = Dict("zigzag" => Int32(0), "up" => Int32(1), "down" => Int32(2),
                         "center" => Int32(3))

function _gen_opts(; counts=(8, 8, 8), size=(1.0, 1.0, 1.0), origin=(0.0, 0.0, 0.0),
                   radius=1.0, height=1.0, radius2=0.5, angle=2pi, max_cells=20_000_000,
                   cell_type="", split="zigzag", order=1, num_twists=1, complete=false,
                   ratio=1.0, seed=0)
    c = Int64[Int64(v) for v in counts]
    length(c) <= 3 || throw(ArgumentError("meshio++: at most three counts"))
    append!(c, fill(Int64(1), 3 - length(c)))
    sz = Cdouble[Cdouble(v) for v in size]
    append!(sz, fill(1.0, 3 - length(sz)))
    o = Cdouble[Cdouble(v) for v in origin]
    append!(o, fill(0.0, 3 - length(o)))
    (length(sz) == 3 && length(o) == 3) ||
        throw(ArgumentError("meshio++: size and origin take two or three entries"))
    _CGenOpts((c[1], c[2], c[3]), (sz[1], sz[2], sz[3]), (o[1], o[2], o[3]),
              Cdouble(radius), Cdouble(height), Int64(max_cells),
              _code(_GEN_CELLS, cell_type, "cell type"), _code(_GEN_SPLITS, split, "split"),
              Int32(order), Int32(0), Cdouble(radius2), Cdouble(angle), Int32(num_twists),
              Int32(complete ? 1 : 0), Cdouble(ratio), UInt64(seed), (Int64(0),))
end

function _generate(sym::Symbol; kwargs...)
    opts = _gen_opts(; kwargs...)
    ptr = ccall(_sym(sym), Ptr{Cvoid}, (Ref{_CGenOpts},), opts)
    Mesh(_check_ptr(ptr))
end

"""
    rectangle(; dims=(1, 1), size=(1.0, 1.0), origin=(0.0, 0.0, 0.0), cell_type="quad",
              split="zigzag", order=1, max_cells=20_000_000, complete=false) -> Mesh

A regular 2-D grid of quads or triangles at `z = origin[3]`, points x fastest. `split`
(`"zigzag"`, `"up"`, `"down"`, `"center"`) lays out the triangles; `order=2` gives
`quad8` / `triangle6` cells. Side regions `xmin`, `xmax`, `ymin`, `ymax`.
"""
rectangle(; dims=(1, 1), size=(1.0, 1.0), origin=(0.0, 0.0, 0.0), cell_type="quad",
          split="zigzag", order=1, max_cells=20_000_000, complete=false) =
    _generate(:mio_rectangle; counts=dims, size=size, origin=origin, cell_type=cell_type,
              split=split, order=order, max_cells=max_cells, complete=complete)

"""
    box(; dims=(1, 1, 1), size=(1.0, 1.0, 1.0), origin=(0.0, 0.0, 0.0),
        cell_type="hexahedron", order=1, max_cells=20_000_000, complete=false) -> Mesh

A regular 3-D grid of hexahedra, tetrahedra (six Kuhn tetrahedra per cell) or wedges.
Side regions `xmin` ... `zmax`.
"""
box(; dims=(1, 1, 1), size=(1.0, 1.0, 1.0), origin=(0.0, 0.0, 0.0), cell_type="hexahedron",
    order=1, max_cells=20_000_000, complete=false) =
    _generate(:mio_box; counts=dims, size=size, origin=origin, cell_type=cell_type,
              order=order, max_cells=max_cells, complete=complete)

"""
    disk(; radius=1.0, n=8, center=(0.0, 0.0, 0.0), cell_type="quad", order=1,
         max_cells=20_000_000, complete=false) -> Mesh

A disk of `n x n` quads (squircle map) or `n` rings of triangles. Side region `rim`.
"""
disk(; radius=1.0, n=8, center=(0.0, 0.0, 0.0), cell_type="quad", order=1,
     max_cells=20_000_000, complete=false) =
    _generate(:mio_disk; radius=radius, counts=(n,), origin=center, cell_type=cell_type,
              order=order, max_cells=max_cells, complete=complete)

"""
    uv_sphere(; radius=1.0, nu=16, nv=8, center=(0.0, 0.0, 0.0), cell_type="triangle",
              order=1, max_cells=20_000_000, complete=false) -> Mesh

A latitude / longitude sphere surface, oriented outward; `cell_type="quad"` gives quads
between the rings and triangle caps.
"""
uv_sphere(; radius=1.0, nu=16, nv=8, center=(0.0, 0.0, 0.0), cell_type="triangle", order=1,
          max_cells=20_000_000, complete=false) =
    _generate(:mio_uv_sphere; radius=radius, counts=(nu, nv), origin=center,
              cell_type=cell_type, order=order, max_cells=max_cells, complete=complete)

"""
    icosa_sphere(; radius=1.0, subdivisions=2, center=(0.0, 0.0, 0.0), order=1,
                 max_cells=20_000_000, complete=false) -> Mesh

A sphere surface of `20 * 4^subdivisions` triangles, oriented outward.
"""
icosa_sphere(; radius=1.0, subdivisions=2, center=(0.0, 0.0, 0.0), order=1,
             max_cells=20_000_000, complete=false) =
    _generate(:mio_icosa_sphere; radius=radius, counts=(subdivisions,), origin=center,
              order=order, max_cells=max_cells, complete=complete)

"""
    cylinder(; radius=1.0, height=1.0, nxy=8, nz=8, center=(0.0, 0.0, 0.0),
             cell_type="hexahedron", order=1, max_cells=20_000_000, complete=false) -> Mesh

A solid cylinder along z, hexahedra or wedges. Side regions `side`, `top`, `bottom`.
"""
cylinder(; radius=1.0, height=1.0, nxy=8, nz=8, center=(0.0, 0.0, 0.0),
         cell_type="hexahedron", order=1, max_cells=20_000_000, complete=false) =
    _generate(:mio_cylinder; radius=radius, height=height, counts=(nxy, nz), origin=center,
              cell_type=cell_type, order=order, max_cells=max_cells, complete=complete)

"""
    reference_triangle(; n=4, size=(1.0, 1.0), origin=(0.0, 0.0, 0.0), order=1,
                       max_cells=20_000_000, complete=false) -> Mesh

The reference triangle `(0,0) (1,0) (0,1)`, scaled by `size` and refined into `n^2`
triangles in barycentric lattice order. Side regions `ymin`, `xmin`, `hypotenuse`.
"""
reference_triangle(; n=4, size=(1.0, 1.0), origin=(0.0, 0.0, 0.0), order=1,
                   max_cells=20_000_000, complete=false) =
    _generate(:mio_reference_triangle; counts=(n,), size=size, origin=origin, order=order,
              max_cells=max_cells, complete=complete)

"""
    ngon(; sides=6, radius=1.0, rings=1, center=(0.0, 0.0, 0.0), order=1,
         max_cells=20_000_000, complete=false) -> Mesh

A regular polygon fanned from its centre (`rings=1`), or refined into `rings` rings of
triangles. Side region `rim`.
"""
ngon(; sides=6, radius=1.0, rings=1, center=(0.0, 0.0, 0.0), order=1,
     max_cells=20_000_000, complete=false) =
    _generate(:mio_ngon; radius=radius, counts=(sides, rings), origin=center, order=order,
              max_cells=max_cells, complete=complete)

"""
    annulus(; inner_radius=0.5, outer_radius=1.0, ntheta=16, nr=4, angle=2pi,
            center=(0.0, 0.0, 0.0), cell_type="quad", order=1,
            max_cells=20_000_000, complete=false) -> Mesh

An annulus, or a sector of `angle` radians whose two end rays stay separate node chains
(the seam pair `match_periodic_nodes` links). Side regions `inner`, `outer` and, for a
sector, `start` and `end`.
"""
annulus(; inner_radius=0.5, outer_radius=1.0, ntheta=16, nr=4, angle=2pi,
        center=(0.0, 0.0, 0.0), cell_type="quad", order=1, max_cells=20_000_000, complete=false) =
    _generate(:mio_annulus; radius=outer_radius, radius2=inner_radius,
              counts=(ntheta, nr), angle=angle, origin=center, cell_type=cell_type,
              order=order, max_cells=max_cells, complete=complete)

"""
    lshape(; n=4, size=2.0, center=(0.0, 0.0, 0.0), cell_type="quad", order=1,
           max_cells=20_000_000, complete=false) -> Mesh

The square of side `size` minus its fourth quadrant (`x > 0`, `y < 0`), the re-entrant
corner at `center`. Side regions `xmin`, `xmax`, `ymin`, `ymax`, `notch_x`, `notch_y`.
"""
lshape(; n=4, size=2.0, center=(0.0, 0.0, 0.0), cell_type="quad", order=1,
       max_cells=20_000_000, complete=false) =
    _generate(:mio_lshape; counts=(n,), size=(size, size, size), origin=center,
              cell_type=cell_type, order=order, max_cells=max_cells, complete=complete)

"""
    plate_with_hole(; radius=0.25, size=1.0, nt=8, nr=4, center=(0.0, 0.0, 0.0),
                    cell_type="quad", order=1, max_cells=20_000_000, complete=false) -> Mesh

A square plate of side `size` around a circular hole, in four quad blocks of `nt x nr`
cells. Side regions `hole`, `xmin`, `xmax`, `ymin`, `ymax`.
"""
plate_with_hole(; radius=0.25, size=1.0, nt=8, nr=4, center=(0.0, 0.0, 0.0),
                cell_type="quad", order=1, max_cells=20_000_000, complete=false) =
    _generate(:mio_plate_with_hole; radius=radius, size=(size, size, size),
              counts=(nt, nr), origin=center, cell_type=cell_type, order=order,
              max_cells=max_cells, complete=complete)

"""
    honeycomb(; dims=(4, 4), edge=1.0, origin=(0.0, 0.0, 0.0), max_cells=20_000_000) -> Mesh

A tiling of pointy-top hexagons as one ragged `polygon` block. Side region `boundary`.
"""
honeycomb(; dims=(4, 4), edge=1.0, origin=(0.0, 0.0, 0.0), max_cells=20_000_000) =
    _generate(:mio_honeycomb; counts=dims, size=(edge, edge, edge), origin=origin,
              max_cells=max_cells)

# --- curves and point clouds ---

"""
    line(; n=8, start=(0.0, 0.0, 0.0), stop=(1.0, 0.0, 0.0), order=1,
         max_cells=20_000_000, complete=false) -> Mesh

A straight segment of `n` `line` cells from `start` to `stop`. Point regions `start`
and `end`.
"""
line(; n=8, start=(0.0, 0.0, 0.0), stop=(1.0, 0.0, 0.0), order=1, max_cells=20_000_000,
     complete=false) =
    _generate(:mio_line; counts=(n,), origin=start, size=stop, order=order,
              max_cells=max_cells, complete=complete)

"""
    circle(; radius=1.0, n=16, center=(0.0, 0.0, 0.0), order=1, max_cells=20_000_000,
           complete=false) -> Mesh

A closed polyline of `n` `line` cells on a circle in the `z = center[3]` plane.
"""
circle(; radius=1.0, n=16, center=(0.0, 0.0, 0.0), order=1, max_cells=20_000_000,
       complete=false) =
    _generate(:mio_circle; radius=radius, counts=(n,), origin=center, order=order,
              max_cells=max_cells, complete=complete)

"""
    helix(; radius=1.0, pitch=1.0, turns=2.0, n=32, center=(0.0, 0.0, 0.0), order=1,
          max_cells=20_000_000, complete=false) -> Mesh

A helix of `n` `line` cells around the z axis, rising `pitch` per turn. Point regions
`start` and `end`.
"""
helix(; radius=1.0, pitch=1.0, turns=2.0, n=32, center=(0.0, 0.0, 0.0), order=1,
      max_cells=20_000_000, complete=false) =
    _generate(:mio_helix; radius=radius, height=pitch, angle=2pi * turns, counts=(n,),
              origin=center, order=order, max_cells=max_cells, complete=complete)

"""
    sphere_points(; n=100, radius=1.0, center=(0.0, 0.0, 0.0), max_cells=20_000_000) -> Mesh

`n` points on a sphere on a Fibonacci lattice, as `vertex` cells.
"""
sphere_points(; n=100, radius=1.0, center=(0.0, 0.0, 0.0), max_cells=20_000_000) =
    _generate(:mio_sphere_points; counts=(n,), radius=radius, origin=center,
              max_cells=max_cells)

"""
    box_points(; n=100, size=(1.0, 1.0, 1.0), origin=(0.0, 0.0, 0.0), seed=0,
               max_cells=20_000_000) -> Mesh

`n` points drawn uniformly in a box from splitmix64 seeded with `seed`, as `vertex`
cells: the same cloud on every platform and binding.
"""
box_points(; n=100, size=(1.0, 1.0, 1.0), origin=(0.0, 0.0, 0.0), seed=0,
           max_cells=20_000_000) =
    _generate(:mio_box_points; counts=(n,), size=size, origin=origin, seed=seed,
              max_cells=max_cells)

# --- surfaces ---

"""
    tetra_sphere(; radius=1.0, subdivisions=3, center=(0.0, 0.0, 0.0), order=1,
                 max_cells=20_000_000, complete=false) -> Mesh

A sphere surface from a regular tetrahedron split in four `subdivisions` times.
"""
tetra_sphere(; radius=1.0, subdivisions=3, center=(0.0, 0.0, 0.0), order=1,
             max_cells=20_000_000, complete=false) =
    _generate(:mio_tetra_sphere; radius=radius, counts=(subdivisions,), origin=center,
              order=order, max_cells=max_cells, complete=complete)

"""
    octa_sphere(; radius=1.0, subdivisions=2, center=(0.0, 0.0, 0.0), order=1,
                max_cells=20_000_000, complete=false) -> Mesh

A sphere surface from a regular octahedron split in four `subdivisions` times.
"""
octa_sphere(; radius=1.0, subdivisions=2, center=(0.0, 0.0, 0.0), order=1,
            max_cells=20_000_000, complete=false) =
    _generate(:mio_octa_sphere; radius=radius, counts=(subdivisions,), origin=center,
              order=order, max_cells=max_cells, complete=complete)

"""
    cube_sphere(; radius=1.0, n=4, center=(0.0, 0.0, 0.0), cell_type="quad", order=1,
                max_cells=20_000_000, complete=false) -> Mesh

A sphere surface from a cube's six faces of `n x n` cells pushed onto the sphere.
"""
cube_sphere(; radius=1.0, n=4, center=(0.0, 0.0, 0.0), cell_type="quad", order=1,
            max_cells=20_000_000, complete=false) =
    _generate(:mio_cube_sphere; radius=radius, counts=(n,), origin=center,
              cell_type=cell_type, order=order, max_cells=max_cells, complete=complete)

"""
    geo_sphere(; radius=1.0, n=4, center=(0.0, 0.0, 0.0), order=1, max_cells=20_000_000,
               complete=false) -> Mesh

A geodesic sphere of frequency `n` (`20 n^2` triangles).
"""
geo_sphere(; radius=1.0, n=4, center=(0.0, 0.0, 0.0), order=1, max_cells=20_000_000,
           complete=false) =
    _generate(:mio_geo_sphere; radius=radius, counts=(n,), origin=center, order=order,
              max_cells=max_cells, complete=complete)

"""
    torus(; major_radius=1.0, minor_radius=0.25, nu=32, nv=16, center=(0.0, 0.0, 0.0),
          cell_type="quad", order=1, max_cells=20_000_000, complete=false) -> Mesh

A torus around the z axis, periodic both ways and welded.
"""
torus(; major_radius=1.0, minor_radius=0.25, nu=32, nv=16, center=(0.0, 0.0, 0.0),
      cell_type="quad", order=1, max_cells=20_000_000, complete=false) =
    _generate(:mio_torus; radius=major_radius, radius2=minor_radius, counts=(nu, nv),
              origin=center, cell_type=cell_type, order=order, max_cells=max_cells,
              complete=complete)

"""
    moebius(; radius=1.0, width=0.5, nu=32, nv=4, num_twists=1, center=(0.0, 0.0, 0.0),
            cell_type="quad", order=1, max_cells=20_000_000, complete=false) -> Mesh

A Moebius band with `num_twists` half-twists (odd: non-orientable). Side region
`boundary`.
"""
moebius(; radius=1.0, width=0.5, nu=32, nv=4, num_twists=1, center=(0.0, 0.0, 0.0),
        cell_type="quad", order=1, max_cells=20_000_000, complete=false) =
    _generate(:mio_moebius; radius=radius, radius2=width, counts=(nu, nv),
              num_twists=num_twists, origin=center, cell_type=cell_type, order=order,
              max_cells=max_cells, complete=complete)

"""
    tube(; radius=1.0, height=1.0, ntheta=16, nz=4, center=(0.0, 0.0, 0.0),
         cell_type="quad", order=1, max_cells=20_000_000, complete=false) -> Mesh

An open cylinder surface around the z axis. Side regions `bottom` and `top`.
"""
tube(; radius=1.0, height=1.0, ntheta=16, nz=4, center=(0.0, 0.0, 0.0), cell_type="quad",
     order=1, max_cells=20_000_000, complete=false) =
    _generate(:mio_tube; radius=radius, height=height, counts=(ntheta, nz), origin=center,
              cell_type=cell_type, order=order, max_cells=max_cells, complete=complete)

"""
    cone(; radius=1.0, height=1.0, ntheta=16, nz=4, center=(0.0, 0.0, 0.0),
         cell_type="quad", order=1, max_cells=20_000_000, complete=false) -> Mesh

A closed cone surface with a sharp rim. Cell regions `side` and `base`.
"""
cone(; radius=1.0, height=1.0, ntheta=16, nz=4, center=(0.0, 0.0, 0.0), cell_type="quad",
     order=1, max_cells=20_000_000, complete=false) =
    _generate(:mio_cone; radius=radius, height=height, counts=(ntheta, nz), origin=center,
              cell_type=cell_type, order=order, max_cells=max_cells, complete=complete)

"""
    capsule(; radius=0.5, height=1.0, ntheta=16, nz=4, ncap=4, center=(0.0, 0.0, 0.0),
            cell_type="quad", order=1, max_cells=20_000_000, complete=false) -> Mesh

A cylinder of length `height` closed by two hemispheres.
"""
capsule(; radius=0.5, height=1.0, ntheta=16, nz=4, ncap=4, center=(0.0, 0.0, 0.0),
        cell_type="quad", order=1, max_cells=20_000_000, complete=false) =
    _generate(:mio_capsule; radius=radius, height=height, counts=(ntheta, nz, ncap),
              origin=center, cell_type=cell_type, order=order, max_cells=max_cells,
              complete=complete)

"""
    box_surface(; dims=(1, 1, 1), size=(1.0, 1.0, 1.0), origin=(0.0, 0.0, 0.0),
                cell_type="quad", order=1, max_cells=20_000_000, complete=false) -> Mesh

The boundary of a box as a closed surface. Cell regions `xmin` ... `zmax`.
"""
box_surface(; dims=(1, 1, 1), size=(1.0, 1.0, 1.0), origin=(0.0, 0.0, 0.0),
            cell_type="quad", order=1, max_cells=20_000_000, complete=false) =
    _generate(:mio_box_surface; counts=dims, size=size, origin=origin, cell_type=cell_type,
              order=order, max_cells=max_cells, complete=complete)

const _PLATONIC_FACES = Dict("tetrahedron" => 4, "cube" => 6, "octahedron" => 8,
                             "dodecahedron" => 12, "icosahedron" => 20)

"""
    platonic(name="icosahedron"; radius=1.0, center=(0.0, 0.0, 0.0),
             max_cells=20_000_000) -> Mesh

A platonic solid's surface, one cell per face: `"tetrahedron"`, `"cube"`,
`"octahedron"`, `"dodecahedron"` (a `polygon` block) or `"icosahedron"`.
"""
function platonic(name::AbstractString="icosahedron"; radius=1.0, center=(0.0, 0.0, 0.0),
                  max_cells=20_000_000)
    haskey(_PLATONIC_FACES, name) ||
        throw(ArgumentError("meshio++: platonic: unknown name '$name'"))
    _generate(:mio_platonic; counts=(_PLATONIC_FACES[name],), radius=radius, origin=center,
              max_cells=max_cells)
end

# --- solids ---

"""
    ball(; radius=1.0, n=4, center=(0.0, 0.0, 0.0), cell_type="hexahedron", order=1,
         max_cells=20_000_000, complete=false) -> Mesh

A ball by the seven-block butterfly, `n` cells from the centre to the surface. Side
region `surface`.
"""
ball(; radius=1.0, n=4, center=(0.0, 0.0, 0.0), cell_type="hexahedron", order=1,
     max_cells=20_000_000, complete=false) =
    _generate(:mio_ball; radius=radius, counts=(n,), origin=center, cell_type=cell_type,
              order=order, max_cells=max_cells, complete=complete)

"""
    spherical_shell(; inner_radius=0.5, outer_radius=1.0, n=4, nr=2, ratio=1.0,
                    center=(0.0, 0.0, 0.0), order=1, max_cells=20_000_000,
                    complete=false) -> Mesh

A thick spherical shell of hexahedra, `ratio` grading the layers. Side regions `inner`
and `outer`.
"""
spherical_shell(; inner_radius=0.5, outer_radius=1.0, n=4, nr=2, ratio=1.0,
                center=(0.0, 0.0, 0.0), order=1, max_cells=20_000_000, complete=false) =
    _generate(:mio_spherical_shell; radius=outer_radius, radius2=inner_radius,
              counts=(n, nr), ratio=ratio, origin=center, order=order, max_cells=max_cells,
              complete=complete)

"""
    pipe(; inner_radius=0.5, outer_radius=1.0, height=1.0, ntheta=16, nr=2, nz=4,
         ratio=1.0, center=(0.0, 0.0, 0.0), order=1, max_cells=20_000_000,
         complete=false) -> Mesh

A thick-walled pipe of hexahedra. Side regions `inner`, `outer`, `bottom`, `top`.
"""
pipe(; inner_radius=0.5, outer_radius=1.0, height=1.0, ntheta=16, nr=2, nz=4, ratio=1.0,
     center=(0.0, 0.0, 0.0), order=1, max_cells=20_000_000, complete=false) =
    _generate(:mio_pipe; radius=outer_radius, radius2=inner_radius, height=height,
              counts=(ntheta, nr, nz), ratio=ratio, origin=center, order=order,
              max_cells=max_cells, complete=complete)

"""
    fichera(; n=2, size=2.0, center=(0.0, 0.0, 0.0), cell_type="hexahedron", order=1,
            max_cells=20_000_000, complete=false) -> Mesh

The cube of side `size` minus its positive octant. Side regions `xmin` ... `zmax`,
`notch_x`, `notch_y`, `notch_z`.
"""
fichera(; n=2, size=2.0, center=(0.0, 0.0, 0.0), cell_type="hexahedron", order=1,
        max_cells=20_000_000, complete=false) =
    _generate(:mio_fichera; counts=(n,), size=(size, size, size), origin=center,
              cell_type=cell_type, order=order, max_cells=max_cells, complete=complete)

"""
    solid_torus(; major_radius=1.0, minor_radius=0.25, n=4, nu=32, center=(0.0, 0.0, 0.0),
                order=1, max_cells=20_000_000, complete=false) -> Mesh

A solid torus of hexahedra, the squircle disk swept around the z axis. Side region
`surface`.
"""
solid_torus(; major_radius=1.0, minor_radius=0.25, n=4, nu=32, center=(0.0, 0.0, 0.0),
            order=1, max_cells=20_000_000, complete=false) =
    _generate(:mio_solid_torus; radius=major_radius, radius2=minor_radius, counts=(n, nu),
              origin=center, order=order, max_cells=max_cells, complete=complete)

"""
    reference_cells(types=String[]; pairs=false) -> Mesh

One cell of each requested meshio type in its reference coordinates, side by side
(every supported type when `types` is empty). `pairs=true` gives two cells sharing a
facet, for linear types.
"""
function reference_cells(types=String[]; pairs::Bool=false)
    names = join(String.(types), ",")
    ptr = ccall(_sym(:mio_reference_cells), Ptr{Cvoid}, (Cstring, Int32), names,
                Int32(pairs))
    Mesh(_check_ptr(ptr))
end

"""
    voxelize(m; resolution=nothing, cell_size=nothing, bounds=nothing, padding=0.0,
             padding_relative=0.0, fill=:all, sign=:pseudonormal,
             watertight_check=:warn, attach_occupancy=false,
             max_cells=20_000_000, plane=:auto) -> NamedTuple

Build a regular grid around `m`. Exactly one of `resolution` and `cell_size`
must be given. `fill` is `:all`, `:surface` or `:inside`. A zero third
`resolution` entry selects 2-D pixels around a planar skin, whose plane is
`plane` (`:auto`, `:xy`, `:xz` or `:yz`; `:auto` takes the axis of zero extent).

Returns `(; mesh, dims, origin, spacing, num_occupied)`.
"""
function voxelize(m::Mesh; resolution=nothing, cell_size=nothing, bounds=nothing,
                  padding::Real=0.0, padding_relative::Real=0.0, fill=:all,
                  sign=:pseudonormal, watertight_check=:warn,
                  attach_occupancy::Bool=false, max_cells::Integer=20_000_000,
                  plane=:auto)
    res = resolution === nothing ? Int64[] : Int64[Int64(v) for v in resolution]
    (resolution === nothing || length(res) == 3) ||
        throw(ArgumentError("resolution must have three entries"))
    bnd = bounds === nothing ? Cdouble[] : Cdouble[Cdouble(v) for v in bounds]
    (bounds === nothing || length(bnd) == 6) ||
        throw(ArgumentError("bounds must have six entries"))

    dims = Vector{Int64}(undef, 3)
    origin = Vector{Cdouble}(undef, 3)
    spacing = Vector{Cdouble}(undef, 3)
    occupied = Ref{Int64}(0)
    ptr = GC.@preserve res bnd dims origin spacing begin
        opts = _CVoxelOpts(isempty(res) ? Ptr{Int64}(C_NULL) : pointer(res),
                           isempty(bnd) ? Ptr{Cdouble}(C_NULL) : pointer(bnd),
                           Cdouble(cell_size === nothing ? 0.0 : cell_size),
                           Cdouble(padding), Cdouble(padding_relative),
                           Int64(max_cells), _code(_VOXEL_FILLS, fill, "fill"),
                           attach_occupancy ? Int32(1) : Int32(0),
                           _code(_SDF_SIGNS, sign, "sign"),
                           _code(_SDF_CHECKS, watertight_check, "watertight check"),
                           _code(_PLANAR_PLANES, plane, "plane"), Int32(0),
                           (Int64(0), Int64(0), Int64(0), Int64(0), Int64(0)))
        ccall(_sym(:mio_voxelize), Ptr{Cvoid},
              (Ptr{Cvoid}, Ref{_CVoxelOpts}, Ptr{Int64}, Ptr{Cdouble}, Ptr{Cdouble},
               Ptr{Int64}),
              _handle(m), Ref(opts), pointer(dims), pointer(origin), pointer(spacing),
              occupied)
    end
    (mesh=Mesh(_check_ptr(ptr)), dims=Tuple(Int.(dims)), origin=Tuple(Float64.(origin)),
     spacing=Tuple(Float64.(spacing)), num_occupied=Int(occupied[]))
end

const _OCTREE_BALANCES = Dict("none" => Int32(0), "face" => Int32(1), "full" => Int32(2))
const _OCTREE_TRANSITIONS = Dict("hanging" => Int32(0), "polyhedra" => Int32(1),
                                 "templates" => Int32(2))
const _PLANAR_PLANES = Dict("auto" => Int32(0), "xy" => Int32(1), "xz" => Int32(2),
                            "yz" => Int32(3))

"""
    octree(m; root_cells=(1, 1, 1), min_level=0, max_level=4, bounds=nothing,
           padding=0.0, padding_relative=0.0, dimension=0, sample_distance=0.0,
           balance=:face, fill=:all, sign=:pseudonormal, watertight_check=:warn,
           attach_distance=false, even_odd=false, max_cells=20_000_000,
           transition=:hanging, snap=false, surface_level=-1, curvature_threshold=0.0,
           feature_angle=0.0, size_array="", plane=:auto) -> NamedTuple

Build an octree (a quadtree for a planar `line` skin) around `m`: leaves are
fine only where the skin is. `balance` is `:none`, `:face` or `:full`; `fill`
is `:all`, `:surface` or `:inside` as for [`voxelize`](@ref). The mesh carries
`octree:level` and `octree:cut` cell data and `octree:hanging` point data.

`transition` is `:hanging`, `:polyhedra` or `:templates` (conforming output; see
`doc/octree.md`). `snap=true` moves the outer nodes onto the skin without
inverting a leaf. `surface_level`, `curvature_threshold`, `feature_angle` and
`size_array` choose how far cut leaves refine. `plane` is `:auto`, `:xy`, `:xz`
or `:yz` for a 2-D skin.

Returns `(; mesh, dimension, num_leaves, num_hanging, finest_level, origin,
root_size, num_snapped, num_refused)`.
"""
function octree(m::Mesh; root_cells=(1, 1, 1), min_level::Integer=0,
                max_level::Integer=4, bounds=nothing, padding::Real=0.0,
                padding_relative::Real=0.0, dimension::Integer=0,
                sample_distance::Real=0.0, balance=:face, fill=:all,
                sign=:pseudonormal, watertight_check=:warn,
                attach_distance::Bool=false, even_odd::Bool=false,
                max_cells::Integer=20_000_000, transition=:hanging, snap::Bool=false,
                surface_level::Integer=-1, curvature_threshold::Real=0.0,
                feature_angle::Real=0.0, size_array::AbstractString="", plane=:auto)
    roots = Int64[Int64(v) for v in root_cells]
    length(roots) == 3 || throw(ArgumentError("root_cells must have three entries"))
    bnd = bounds === nothing ? Cdouble[] : Cdouble[Cdouble(v) for v in bounds]
    (bounds === nothing || length(bnd) == 6) ||
        throw(ArgumentError("bounds must have six entries"))

    report = Ref{_COctreeReport}()
    size_name = String(size_array)
    ptr = GC.@preserve roots bnd size_name begin
        opts = _COctreeOpts(pointer(roots),
                            isempty(bnd) ? Ptr{Cdouble}(C_NULL) : pointer(bnd),
                            Cdouble(padding), Cdouble(padding_relative),
                            Cdouble(sample_distance), Int64(min_level), Int64(max_level),
                            Int64(dimension), Int64(max_cells),
                            _code(_OCTREE_BALANCES, balance, "balance"),
                            _code(_VOXEL_FILLS, fill, "fill"),
                            attach_distance ? Int32(1) : Int32(0),
                            even_odd ? Int32(1) : Int32(0),
                            _code(_SDF_SIGNS, sign, "sign"),
                            _code(_SDF_CHECKS, watertight_check, "watertight check"),
                            isempty(size_name) ? Cstring(C_NULL) : Base.unsafe_convert(Cstring, size_name),
                            Cdouble(curvature_threshold), Cdouble(feature_angle),
                            Int64(surface_level),
                            _code(_OCTREE_TRANSITIONS, transition, "transition"),
                            snap ? Int32(1) : Int32(0),
                            _code(_PLANAR_PLANES, plane, "plane"), Int32(0),
                            (Int64(0), Int64(0), Int64(0), Int64(0)))
        ccall(_sym(:mio_octree), Ptr{Cvoid},
              (Ptr{Cvoid}, Ref{_COctreeOpts}, Ref{_COctreeReport}),
              _handle(m), Ref(opts), report)
    end
    mesh = Mesh(_check_ptr(ptr))
    r = report[]
    (mesh=mesh, dimension=Int(r.dimension), num_leaves=Int(r.num_leaves),
     num_hanging=Int(r.num_hanging), finest_level=Int(r.finest_level),
     origin=Tuple(Float64.(r.origin)), root_size=Float64(r.root_size),
     num_snapped=Int(r.num_snapped), num_refused=Int(r.num_refused))
end

const _OCTREE_MESH_OUTPUTS = Dict("castellated" => Int32(0), "tet" => Int32(1),
                                  "hex" => Int32(2))
const _OCTREE_MESH_OUTPUT_NAMES = ("castellated", "tet", "hex")
const _OCTREE_MESH_CONFORMITIES = Dict("stuffing" => Int32(0), "cleaving" => Int32(1))
const _OCTREE_MESH_SKIN_GRADINGS = Dict("graded" => Int32(0), "uniform" => Int32(1))
const _OCTREE_MESH_FILLS = Dict("all" => Int32(0), "inside" => Int32(1), "cover" => Int32(2))
const _OCTREE_MESH_HEX_GRADINGS = Dict("uniform" => Int32(0), "graded" => Int32(1))

"""
    octree_mesh(m; output=:tet, min_size=0.0, max_size=0.0, max_level=-1, bounds=nothing,
                padding=0.0, padding_relative=0.1, curvature_angle=15.0,
                feature_angle=30.0, corner_angle=30.0, feature_size=0.0, gap_cells=2.0,
                gap_outside=false, size_array=nothing, surface_region=nothing,
                grading_layers=1, skin_grading=:graded, surface_size=0.0,
                winding_threshold=0.5, winding_beta=2.0, conformity=:stuffing,
                alpha_long=0.24799, alpha_short=0.40173, preserve_features=true,
                fit_iterations=5, dihedral_floor_guard=true, fill=:inside,
                keep_outside=false, attach_distance=true, attach_quality=false,
                pillow=true, hex_grading=:uniform, untangle_iterations=50,
                max_cells=20_000_000, max_tets=100_000_000, volume_fix=true,
                zero_volume_factor=1e-6, cap_flatness=0.01, crossing_margin=0.01,
                crossing_gap=0.05) -> NamedTuple

Mesh the inside of the skin `m` on a balanced octree: conforming graded
tetrahedra by default (`output=:tet`), or the castellated leaves
(`output=:castellated`) with surrogate-boundary regions and `sdf:distance` /
`sdf:closest_point`; `:hex` is the all-hexahedral dual. `conformity` is
`:stuffing` or `:cleaving`, `skin_grading` `:graded` or `:uniform`, `fill`
(castellated output) `:all`, `:inside` or `:cover`, and `hex_grading`
(hexahedral output) `:uniform` or `:graded`; an unknown one throws
`ArgumentError`. Every keyword left at `nothing` takes its default from
`mio_octree_mesh_opts_init`, i.e. from the C side; the defaults shown above are
those. See `doc/octree_mesh.md` for every option and report field.

Returns `(; mesh, output, guarantee_applies, origin, root_size, root_cells,
max_level, finest_level, surface_level, min_size, max_size, num_leaves,
num_points, num_cells, num_hexahedra, num_tetrahedra, num_inside,
num_boundary_cells, num_surrogate_faces, num_hanging, num_cut_points,
num_warped, num_warp_refused, min_dihedral, max_dihedral, min_scaled_jacobian,
hausdorff_to_skin, num_inverted, num_boundary_edges_open,
num_non_manifold_boundary_edges, num_corners, num_corners_snapped,
num_feature_chains, num_feature_chains_missed, num_feature_moves,
num_optimized_moves, num_flips, num_winding_rechecks, num_winding_skipped, num_gaps_filled, num_gaps_skipped, num_collapsed, num_reordered, num_overlaps, num_caps_removed, num_split,
num_slivers, num_untangled, num_warnings, leaves_per_level, timings, warnings)`,
with `output` the name of the output produced and `warnings` a
`Vector{String}`, one entry per warning.
"""
function octree_mesh(m::Mesh; output=nothing, min_size=nothing, max_size=nothing,
                     max_level=nothing, bounds=nothing, padding=nothing,
                     padding_relative=nothing, curvature_angle=nothing,
                     feature_angle=nothing, corner_angle=nothing, feature_size=nothing,
                     gap_cells=nothing, gap_outside=nothing, size_array=nothing,
                     surface_region=nothing, grading_layers=nothing, skin_grading=nothing,
                     surface_size=nothing, winding_threshold=nothing, winding_beta=nothing,
                     conformity=nothing, alpha_long=nothing, alpha_short=nothing,
                     preserve_features=nothing, fit_iterations=nothing,
                     dihedral_floor_guard=nothing, fill=nothing, keep_outside=nothing,
                     attach_distance=nothing, attach_quality=nothing, pillow=nothing,
                     hex_grading=nothing, untangle_iterations=nothing, max_cells=nothing,
                     max_tets=nothing, volume_fix=nothing, zero_volume_factor=nothing,
                     cap_flatness=nothing, crossing_margin=nothing, crossing_gap=nothing,
                     post_adapt=nothing, adapt_hmin=nothing, adapt_hmax=nothing,
                     adapt_hgrad=nothing, adapt_hausd=nothing, adapt_max_iterations=nothing,
                     adapt_preserve_boundary=nothing)
    bnd = bounds === nothing ? Cdouble[] : Cdouble[Cdouble(v) for v in bounds]
    (bounds === nothing || length(bnd) == 6) ||
        throw(ArgumentError("bounds must have six entries"))
    size_name = size_array === nothing ? "" : String(size_array)
    region_name = surface_region === nothing ? "" : String(surface_region)

    opts = Ref{_COctreeMeshOpts}()
    ccall(_sym(:mio_octree_mesh_opts_init), Cvoid, (Ref{_COctreeMeshOpts},), opts)
    output === nothing ||
        _setopt!(opts, :output, _code(_OCTREE_MESH_OUTPUTS, output, "output"))
    conformity === nothing ||
        _setopt!(opts, :conformity, _code(_OCTREE_MESH_CONFORMITIES, conformity, "conformity"))
    skin_grading === nothing ||
        _setopt!(opts, :skin_grading,
                 _code(_OCTREE_MESH_SKIN_GRADINGS, skin_grading, "skin grading"))
    fill === nothing || _setopt!(opts, :fill, _code(_OCTREE_MESH_FILLS, fill, "fill"))
    hex_grading === nothing ||
        _setopt!(opts, :hex_grading,
                 _code(_OCTREE_MESH_HEX_GRADINGS, hex_grading, "hex grading"))
    scalars = (; min_size, max_size, padding, padding_relative, curvature_angle,
               feature_angle, corner_angle, feature_size, gap_cells, surface_size,
               winding_threshold, winding_beta, alpha_long, alpha_short, max_level,
               grading_layers, fit_iterations, untangle_iterations, max_cells, max_tets,
               zero_volume_factor, cap_flatness, crossing_margin, crossing_gap)
    for (k, v) in pairs(scalars)
        v === nothing || _setopt!(opts, k, v)
    end
    flags = (; gap_outside, preserve_features, dihedral_floor_guard, keep_outside,
             attach_distance, attach_quality, pillow, volume_fix)
    for (k, v) in pairs(flags)
        v === nothing || _setopt!(opts, k, Bool(v) ? Int32(1) : Int32(0))
    end

    adapt_opts = Ref{_CAdaptVolumeOpts}()
    _adapt_volume_opts!(adapt_opts; hmin=adapt_hmin, hmax=adapt_hmax, hgrad=adapt_hgrad,
                        hausd=adapt_hausd, max_iterations=adapt_max_iterations,
                        preserve_boundary=adapt_preserve_boundary)
    post_adapt === nothing ||
        _setopt!(opts, :post_adapt, Bool(post_adapt) ? Int32(1) : Int32(0))
    report = Ref{_COctreeMeshReport}()
    ptr = GC.@preserve bnd size_name region_name adapt_opts begin
        _setopt!(opts, :adapt,
                 Ptr{Cvoid}(Base.unsafe_convert(Ptr{_CAdaptVolumeOpts}, adapt_opts)))
        isempty(bnd) || _setopt!(opts, :bounds, pointer(bnd))
        isempty(size_name) ||
            _setopt!(opts, :size_array, Base.unsafe_convert(Cstring, size_name))
        isempty(region_name) ||
            _setopt!(opts, :surface_region, Base.unsafe_convert(Cstring, region_name))
        ccall(_sym(:mio_octree_mesh), Ptr{Cvoid},
              (Ptr{Cvoid}, Ref{_COctreeMeshOpts}, Ref{_COctreeMeshReport}),
              _handle(m), opts, report)
    end
    mesh = Mesh(_check_ptr(ptr))
    text = _getstring((buf, n) -> ccall(_sym(:mio_octree_mesh_warnings), Int64,
                                        (Ptr{UInt8}, Int64), buf, n))
    stage_text = _getstring((buf, n) -> ccall(_sym(:mio_octree_mesh_stage_names), Int64,
                                              (Ptr{UInt8}, Int64), buf, n))
    stage_count = ccall(_sym(:mio_octree_mesh_stage_seconds), Int64,
                        (Ptr{Cdouble}, Int64), C_NULL, 0)
    stage_seconds = Vector{Cdouble}(undef, max(stage_count, 0))
    isempty(stage_seconds) ||
        ccall(_sym(:mio_octree_mesh_stage_seconds), Int64, (Ptr{Cdouble}, Int64),
              stage_seconds, length(stage_seconds))
    level_count = ccall(_sym(:mio_octree_mesh_leaves_per_level), Int64,
                        (Ptr{Int64}, Int64), C_NULL, 0)
    leaves_per_level = Vector{Int64}(undef, max(level_count, 0))
    isempty(leaves_per_level) ||
        ccall(_sym(:mio_octree_mesh_leaves_per_level), Int64, (Ptr{Int64}, Int64),
              leaves_per_level, length(leaves_per_level))
    stages = isempty(stage_text) ? String[] : String.(Base.split(stage_text, '\n'))
    adapt_report = Ref{_CAdaptVolumeReport}()
    adapted = ccall(_sym(:mio_octree_mesh_adapt_report), Cint, (Ref{_CAdaptVolumeReport},),
                    adapt_report) == 1
    r = report[]
    (mesh=mesh,
     output=0 <= r.output < length(_OCTREE_MESH_OUTPUT_NAMES) ?
            _OCTREE_MESH_OUTPUT_NAMES[r.output + 1] : string(r.output),
     guarantee_applies=r.guarantee_applies != 0, origin=Tuple(Float64.(r.origin)),
     root_size=Float64(r.root_size), root_cells=Tuple(Int.(r.root_cells)),
     max_level=Int(r.max_level), finest_level=Int(r.finest_level),
     surface_level=Int(r.surface_level), min_size=Float64(r.min_size),
     max_size=Float64(r.max_size), num_leaves=Int(r.num_leaves),
     num_points=Int(r.num_points), num_cells=Int(r.num_cells),
     num_hexahedra=Int(r.num_hexahedra), num_tetrahedra=Int(r.num_tetrahedra),
     num_inside=Int(r.num_inside), num_boundary_cells=Int(r.num_boundary_cells),
     num_surrogate_faces=Int(r.num_surrogate_faces), num_hanging=Int(r.num_hanging),
     num_cut_points=Int(r.num_cut_points), num_warped=Int(r.num_warped),
     num_warp_refused=Int(r.num_warp_refused), min_dihedral=Float64(r.min_dihedral),
     max_dihedral=Float64(r.max_dihedral),
     min_scaled_jacobian=Float64(r.min_scaled_jacobian),
     hausdorff_to_skin=Float64(r.hausdorff_to_skin), num_inverted=Int(r.num_inverted),
     num_boundary_edges_open=Int(r.num_boundary_edges_open),
     num_non_manifold_boundary_edges=Int(r.num_non_manifold_boundary_edges),
     num_corners=Int(r.num_corners), num_corners_snapped=Int(r.num_corners_snapped),
     num_feature_chains=Int(r.num_feature_chains),
     num_feature_chains_missed=Int(r.num_feature_chains_missed),
     num_feature_moves=Int(r.num_feature_moves),
     num_optimized_moves=Int(r.num_optimized_moves), num_flips=Int(r.num_flips),
     num_winding_rechecks=Int(r.num_winding_rechecks),
     num_winding_skipped=Int(r.num_winding_skipped),
     num_gaps_filled=Int(r.num_gaps_filled), num_gaps_skipped=Int(r.num_gaps_skipped),
     num_collapsed=Int(r.num_collapsed), num_reordered=Int(r.num_reordered),
     num_overlaps=Int(r.num_overlaps), num_caps_removed=Int(r.num_caps_removed),
     num_split=Int(r.num_split), num_slivers=Int(r.num_slivers),
     num_untangled=Int(r.num_untangled),
     num_warnings=Int(r.num_warnings),
     leaves_per_level=Int.(leaves_per_level),
     timings=[(stages[i], Float64(stage_seconds[i])) for i in 1:min(length(stages), length(stage_seconds))],
     warnings=isempty(text) ? String[] : String.(Base.split(text, '\n')),
     adapt=adapted ? _adapt_volume_named(adapt_report[]) : nothing)
end

const _CURVATURE_DUAL_AREAS = Dict(:mixed_voronoi => Int32(0), :barycentric => Int32(1))

"""
    compute_curvature(m; mean=true, gaussian=true, dual_area=:mixed_voronoi,
                      include_boundary=false, record_area=false,
                      record_principal=false, region="")
        -> (; mesh, quality, num_boundary, num_isolated, num_degenerate,
             total_angle_defect)

Per-vertex mean (`H`) and Gaussian (`K`) curvature of a surface mesh, by the
angle defect for `K` and the cotangent Laplace-Beltrami operator for `H` — the
signed distance's natural companion as a node feature.

Writes `curvature:mean` and `curvature:gaussian` as point data, optionally
`curvature:area` (the dual area each was divided by) and `curvature:principal`
(`(n, 2)`, `k1 >= k2`). Geometry, connectivity and existing data are carried
through unchanged.

`total_angle_defect` is the oracle: on a CLOSED surface it is `2*pi*chi`
exactly — `4*pi` for anything sphere-like — whatever the tessellation and
whichever `dual_area`, so a value that is not that means the input is not
closed or the result is not sane.

`H` is orientation-dependent and `K` is not, so check
`quality.inconsistent_pairs` before trusting a sign: a nonzero count means
facets disagree about which side is out. This never repairs its input.
"""
function compute_curvature(m::Mesh; mean::Bool=true, gaussian::Bool=true,
                           dual_area::Symbol=:mixed_voronoi,
                           include_boundary::Bool=false, record_area::Bool=false,
                           record_principal::Bool=false, region::AbstractString="")
    haskey(_CURVATURE_DUAL_AREAS, dual_area) ||
        throw(ArgumentError("meshio++: curvature: unknown dual area '$(dual_area)' " *
                            "(expected :mixed_voronoi or :barycentric)"))
    report = Ref{_CCurvatureReport}()
    region_c = Vector{UInt8}(codeunits(String(region) * "\0"))
    ptr = GC.@preserve region_c begin
        opts = _CCurvatureOpts(Cstring(pointer(region_c)),
                               mean ? Int32(1) : Int32(0),
                               gaussian ? Int32(1) : Int32(0),
                               _CURVATURE_DUAL_AREAS[dual_area],
                               include_boundary ? Int32(1) : Int32(0),
                               record_area ? Int32(1) : Int32(0),
                               record_principal ? Int32(1) : Int32(0),
                               (Int64(0), Int64(0), Int64(0), Int64(0), Int64(0), Int64(0)))
        ccall(_sym(:mio_compute_curvature), Ptr{Cvoid},
              (Ptr{Cvoid}, Ref{_CCurvatureOpts}, Ptr{_CCurvatureReport}),
              _handle(m), Ref(opts), report)
    end
    r = _check_ptr(ptr)
    rep = report[]
    q = rep.quality
    (mesh=Mesh(r),
     quality=(boundary_edges=Int(q.boundary_edges),
              non_manifold_edges=Int(q.non_manifold_edges),
              inconsistent_pairs=Int(q.inconsistent_pairs),
              degenerate_triangles=Int(q.degenerate_triangles),
              watertight=q.watertight != 0),
     num_boundary=Int(rep.num_boundary), num_isolated=Int(rep.num_isolated),
     num_degenerate=Int(rep.num_degenerate),
     total_angle_defect=Float64(rep.total_angle_defect))
end

const _NORMALS_WEIGHTS = Dict(:angle => Int32(0), :area => Int32(1))

"""
    compute_normals(m; point_normals=true, cell_normals=false, weight=:angle,
                    split_angle=nothing, record_parent_ids=false, region="")
        -> (; mesh, quality, num_isolated, num_undefined, num_degenerate,
             num_split_points, num_added_points)

Point and/or cell normals of a surface mesh, optionally splitting vertices at
creases so every point carries exactly one normal.

Writes `normals` as point data `(n, 3)` and, with `cell_normals`, as cell data.
`split_angle=nothing` gives one smooth normal per point; a number of degrees in
`[0, 180]` duplicates points wherever the surface creases by more than that,
appending the copies after the original points (cells keep their numbering).
`record_parent_ids` adds `normals:parent_point`. `weight` is `:angle` (default)
or `:area`.

Never reorients: check `quality.inconsistent_pairs` before trusting a normal; a
nonzero count means some normals average faces that disagree about which side is
out. A volume block is refused by name pointing at `extract_surface`, a
higher-order one pointing at `convert_cells`' `:linearize`.
"""
function compute_normals(m::Mesh; point_normals::Bool=true, cell_normals::Bool=false,
                         weight::Symbol=:angle,
                         split_angle::Union{Nothing,Real}=nothing,
                         record_parent_ids::Bool=false, region::AbstractString="")
    haskey(_NORMALS_WEIGHTS, weight) ||
        throw(ArgumentError("meshio++: normals: unknown weight '$(weight)' " *
                            "(expected :angle or :area)"))
    report = Ref{_CNormalsReport}()
    region_c = Vector{UInt8}(codeunits(String(region) * "\0"))
    ptr = GC.@preserve region_c begin
        opts = _CNormalsOpts(Cstring(pointer(region_c)),
                             point_normals ? Int32(1) : Int32(0),
                             cell_normals ? Int32(1) : Int32(0),
                             _NORMALS_WEIGHTS[weight],
                             split_angle === nothing ? Int32(0) : Int32(1),
                             split_angle === nothing ? 30.0 : Float64(split_angle),
                             record_parent_ids ? Int32(1) : Int32(0), Int32(0),
                             (Int64(0), Int64(0), Int64(0), Int64(0), Int64(0), Int64(0)))
        ccall(_sym(:mio_compute_normals), Ptr{Cvoid},
              (Ptr{Cvoid}, Ref{_CNormalsOpts}, Ptr{_CNormalsReport}),
              _handle(m), Ref(opts), report)
    end
    r = _check_ptr(ptr)
    rep = report[]
    (mesh=Mesh(r), quality=_quality_tuple(rep.quality),
     num_isolated=Int(rep.num_isolated), num_undefined=Int(rep.num_undefined),
     num_degenerate=Int(rep.num_degenerate), num_split_points=Int(rep.num_split_points),
     num_added_points=Int(rep.num_added_points))
end

"""
    feature_edges(m; feature_angle=30.0, feature=true, boundary=true,
                  non_manifold=true, inconsistent=true, region="")
        -> (; mesh, num_feature, num_boundary, num_non_manifold, num_inconsistent)

The sharp, open, non-manifold and inconsistently wound edges of a surface (or of
a volume mesh's skin) as a mesh of `line` cells over the input's points, with
cell data `feature:kind` (1 feature, 2 boundary, 3 non-manifold,
4 inconsistent) and `feature:angle` (degrees). The counts are over the whole
surface, whether or not a category was selected. See `doc/feature_edges.md`.
"""
function feature_edges(m::Mesh; feature_angle::Real=30.0, feature::Bool=true,
                       boundary::Bool=true, non_manifold::Bool=true,
                       inconsistent::Bool=true, region::AbstractString="")
    report = Ref{_CFeatureEdgesReport}()
    region_c = Vector{UInt8}(codeunits(String(region) * "\0"))
    ptr = GC.@preserve region_c begin
        opts = _CFeatureEdgesOpts(Cstring(pointer(region_c)), Float64(feature_angle),
                                  Int32(feature), Int32(boundary), Int32(non_manifold),
                                  Int32(inconsistent), ntuple(_ -> Int64(0), 6))
        ccall(_sym(:mio_feature_edges), Ptr{Cvoid},
              (Ptr{Cvoid}, Ref{_CFeatureEdgesOpts}, Ptr{_CFeatureEdgesReport}),
              _handle(m), Ref(opts), report)
    end
    r = _check_ptr(ptr)
    rep = report[]
    (mesh=Mesh(r), num_feature=Int(rep.num_feature), num_boundary=Int(rep.num_boundary),
     num_non_manifold=Int(rep.num_non_manifold), num_inconsistent=Int(rep.num_inconsistent))
end

"""
    segment(m; feature_angle=30.0, normal_tolerance=0.0, curvature_tolerance=0.0,
            min_patch_cells=1, region="", prefix="patch", method="angle",
            tolerance=1e-3)
        -> (; mesh, num_patches, num_merged, num_feature_edges, num_boundary_edges,
             num_non_manifold_edges)

Split a surface into smooth patches by region growing, bounded by its feature
edges. The mesh comes back unchanged plus cell data `segment:patch` (-1 off the
surface), field data `segment:adjacency` (rows `(a, b, shared edges)`) and one
cell region per patch, `<prefix>:<k>`. `method="primitive"` (v17.9.0) ends a
patch where the surface its facets lie on changes, at `tolerance` (relative to
the diagonal), so a fillet comes apart from its tangent planes. See
`doc/segment.md`.
"""
function segment(m::Mesh; feature_angle::Real=30.0, normal_tolerance::Real=0.0,
                 curvature_tolerance::Real=0.0, min_patch_cells::Integer=1,
                 region::AbstractString="", prefix::AbstractString="patch",
                 method::AbstractString="angle", tolerance::Real=1e-3)
    method_c = method == "angle" ? Int32(0) : method == "primitive" ? Int32(1) :
               throw(ArgumentError("meshio++: segment: method must be \"angle\" or \"primitive\""))
    report = Ref{_CSegmentReport}()
    region_c = Vector{UInt8}(codeunits(String(region) * "\0"))
    prefix_c = Vector{UInt8}(codeunits(String(prefix) * "\0"))
    ptr = GC.@preserve region_c prefix_c begin
        opts = _CSegmentOpts(Cstring(pointer(region_c)), Cstring(pointer(prefix_c)),
                             Float64(feature_angle), Float64(normal_tolerance),
                             Float64(curvature_tolerance), Int64(min_patch_cells),
                             method_c, Int32(0), Float64(tolerance), ntuple(_ -> Int64(0), 4))
        ccall(_sym(:mio_segment), Ptr{Cvoid},
              (Ptr{Cvoid}, Ref{_CSegmentOpts}, Ptr{_CSegmentReport}),
              _handle(m), Ref(opts), report)
    end
    r = _check_ptr(ptr)
    rep = report[]
    (mesh=Mesh(r), num_patches=Int(rep.num_patches), num_merged=Int(rep.num_merged),
     num_feature_edges=Int(rep.num_feature_edges),
     num_boundary_edges=Int(rep.num_boundary_edges),
     num_non_manifold_edges=Int(rep.num_non_manifold_edges))
end

"""
    select_region(m, name; kind="cell", box=nothing, plane=nothing, sphere=nothing,
                  where_array="", where_op="<", where_value=0.0, quality_metric="",
                  quality_op="<", quality_value=0.0, component=nothing, seed=nothing,
                  feature_angle=30.0, normal_tolerance=0.0, mode="all",
                  boundary_only=true, tag=-1)
        -> (; mesh, num_candidates, num_selected, num_grown, num_components)

Select points, cells or sides and write them as one named region; the mesh comes
back unchanged plus that region (a region of the same name and kind is
replaced). `kind` is `"cell"`, `"point"` or `"side"`. Every criterion given must
pass: `box = (lo, hi)`, `plane = (point, normal)`, `sphere = (center, radius)`,
`where_array`/`where_op`/`where_value` on point or cell data (as `refine` spells
it), `quality_metric`/`quality_op`/`quality_value` on a `compute_quality` metric, `component = k` (1-based, in the order
`split(by="component")` returns its pieces) and `seed = (x, y, z)`, a patch grown from the nearest facet across
edges no sharper than `feature_angle`. `num_grown` and `num_components` are -1
when unused. See `doc/select.md`.
"""
function select_region(m::Mesh, name::AbstractString; kind::AbstractString="cell",
                       box=nothing, plane=nothing, sphere=nothing,
                       where_array::AbstractString="", where_op::AbstractString="<",
                       where_value::Real=0.0, quality_metric::AbstractString="",
                       quality_op::AbstractString="<", quality_value::Real=0.0,
                       component=nothing, seed=nothing,
                       feature_angle::Real=30.0, normal_tolerance::Real=0.0,
                       mode::AbstractString="all", boundary_only::Bool=true, tag::Integer=-1)
    kinds = ("point", "cell", "side")
    k = findfirst(==(String(kind)), kinds)
    k === nothing && throw(ArgumentError("meshio++: select_region: kind must be \"point\", " *
                                         "\"cell\" or \"side\""))
    modes = ("all", "any", "centroid")
    md = findfirst(==(String(mode)), modes)
    md === nothing && throw(ArgumentError("meshio++: select_region: mode must be \"all\", " *
                                          "\"any\" or \"centroid\""))
    v3(x) = (Float64(x[1]), Float64(x[2]), Float64(x[3]))
    zero3 = (0.0, 0.0, 0.0)
    cmp(op) = get(_REFINE_COMPARES, String(op)) do
        throw(ArgumentError("meshio++: select_region: unknown comparison $(op)"))
    end
    name_c = Vector{UInt8}(codeunits(String(name) * "\0"))
    where_c = Vector{UInt8}(codeunits(String(where_array) * "\0"))
    quality_c = Vector{UInt8}(codeunits(String(quality_metric) * "\0"))
    report = Ref{_CSelectReport}()
    ptr = GC.@preserve name_c where_c quality_c begin
        opts = _CSelectOpts(Cstring(pointer(name_c)), Int32(k - 1), Int32(md - 1), Int64(tag),
                            Int32(boundary_only), Int32(box !== nothing),
                            Int32(plane !== nothing), Int32(sphere !== nothing),
                            Int32(seed !== nothing),
                            cmp(where_op), cmp(quality_op), Int32(0),
                            box === nothing ? zero3 : v3(box[1]),
                            box === nothing ? zero3 : v3(box[2]),
                            plane === nothing ? zero3 : v3(plane[1]),
                            plane === nothing ? (0.0, 0.0, 1.0) : v3(plane[2]),
                            sphere === nothing ? zero3 : v3(sphere[1]),
                            sphere === nothing ? 0.0 : Float64(sphere[2]),
                            Cstring(pointer(where_c)),
                            Float64(where_value),
                            Cstring(pointer(quality_c)),
                            Float64(quality_value),
                            component === nothing ? Int64(-1) : Int64(component) - 1,
                            seed === nothing ? zero3 : v3(seed), Float64(feature_angle),
                            Float64(normal_tolerance), ntuple(_ -> Int64(0), 8))
        ccall(_sym(:mio_select), Ptr{Cvoid}, (Ptr{Cvoid}, Ref{_CSelectOpts}, Ptr{_CSelectReport}),
              _handle(m), Ref(opts), report)
    end
    r = _check_ptr(ptr)
    rep = report[]
    (mesh=Mesh(r), num_candidates=Int(rep.num_candidates), num_selected=Int(rep.num_selected),
     num_grown=Int(rep.num_grown), num_components=Int(rep.num_components))
end

const _FIT_TYPES = ("plane", "sphere", "cylinder", "cone", "torus")

# Every primitive, `(1 << 1) | ... | (1 << 5)` (C `MIO_FIT_ALL_TYPES`).
const _FIT_ALL_TYPES = UInt32(0x3E)

# The default RANSAC seed of `refine_patches`: the bytes of "meshio++".
const _REFINE_DEFAULT_SEED = 0x6d657368696f2b2b

# `mio_refine_stop` codes, in order.
const _REFINE_STOPS = (:unchanged, :small_gain, :iteration_cap, :energy_rose)

# The `1 << type` bit mask of a collection of `_FIT_TYPES` names.
function _fit_type_mask(types, op::AbstractString)
    mask = UInt32(0)
    for t in types
        k = findfirst(==(String(t)), _FIT_TYPES)
        k === nothing && throw(ArgumentError("meshio++: $(op): unknown type $(t)"))
        mask |= UInt32(1) << k
    end
    mask
end

"""
    fit_primitives(m; patch_data="segment:patch", tolerance=1e-3, relative=true,
                   min_inlier_fraction=0.9,
                   types=("plane", "sphere", "cylinder", "cone", "torus"),
                   max_iterations=50, region="", seeds=nothing, feature_angle=30.0,
                   rim_angle=15.0, output_region="fit:feature", used_points="auto",
                   criterion="least_squares")
        -> (; mesh, num_patches, num_plane, num_sphere, num_cylinder, num_cone,
             num_torus, num_unfitted, num_feature_cells, warnings)

Fit a plane, sphere, cylinder, cone or torus to each patch of a surface (the
first accepted in that order, or none). The mesh comes back unchanged plus cell
data `fit:type`, `fit:rms`, `fit:max` and per-patch field data `fit:patch_*`.
`patch_data=""` fits the whole surface as one patch. `region` fits a named cell
region, `seeds` (a vector of 3-vectors) grows and fits a feature from each and
writes its cells as `output_region`, and `used_points` (`"auto"`, `"all"`,
`"3sigma"`, `"2sigma"`, `"1sigma"`) trims outliers (roadmap 7.1.1, 7.1.2);
`criterion` (`"minimum_zone"`, `"min_circumscribed"`, `"max_inscribed"`) refits
planes, spheres and cylinders under that reference (roadmap 7.1.4). See
`doc/fit_primitives.md`.
"""
function fit_primitives(m::Mesh; patch_data::AbstractString="segment:patch",
                        tolerance::Real=1e-3, relative::Bool=true,
                        min_inlier_fraction::Real=0.9,
                        types=("plane", "sphere", "cylinder", "cone", "torus"),
                        max_iterations::Integer=50, region::AbstractString="",
                        seeds=nothing, feature_angle::Real=30.0, rim_angle::Real=15.0,
                        output_region::AbstractString="fit:feature",
                        used_points::AbstractString="auto",
                        criterion::AbstractString="least_squares")
    used = findfirst(==(String(used_points)), ("auto", "all", "3sigma", "2sigma", "1sigma"))
    used === nothing &&
        throw(ArgumentError("meshio++: fit_primitives: unknown used_points $(used_points)"))
    crit = findfirst(==(String(criterion)),
                     ("least_squares", "minimum_zone", "min_circumscribed", "max_inscribed"))
    crit === nothing &&
        throw(ArgumentError("meshio++: fit_primitives: unknown criterion $(criterion)"))
    flat = seeds === nothing ? Cdouble[] : Cdouble[Float64(x) for s in seeds for x in s]
    length(flat) % 3 == 0 ||
        throw(ArgumentError("meshio++: fit_primitives: seeds must be 3-vectors"))
    mask = UInt32(0)
    for t in types
        k = findfirst(==(String(t)), _FIT_TYPES)
        k === nothing && throw(ArgumentError("meshio++: fit_primitives: unknown type $(t)"))
        mask |= UInt32(1) << k
    end
    report = Ref{_CFitPrimitivesReport}()
    patch_c = Vector{UInt8}(codeunits(String(patch_data) * "\0"))
    region_c = Vector{UInt8}(codeunits(String(region) * "\0"))
    out_c = Vector{UInt8}(codeunits(String(output_region) * "\0"))
    ptr = GC.@preserve patch_c region_c out_c flat begin
        base = _CFitPrimitivesOpts(Cstring(pointer(patch_c)), Float64(tolerance),
                                   Float64(min_inlier_fraction), Int32(relative), mask,
                                   Int64(max_iterations), ntuple(_ -> Int64(0), 6))
        opts = _CFitPrimitivesExOpts(base, Cstring(pointer(region_c)),
                                     isempty(flat) ? Ptr{Cdouble}(C_NULL) : pointer(flat),
                                     Int64(length(flat) ÷ 3), Float64(feature_angle),
                                     Float64(rim_angle), Cstring(pointer(out_c)),
                                     Int32(used - 1), Int32(crit - 1), ntuple(_ -> Int64(0), 6))
        ccall(_sym(:mio_fit_primitives_ex), Ptr{Cvoid},
              (Ptr{Cvoid}, Ref{_CFitPrimitivesExOpts}, Ptr{_CFitPrimitivesReport}),
              _handle(m), Ref(opts), report)
    end
    r = _check_ptr(ptr)
    rep = report[]
    text = _getstring((buf, n) -> ccall(_sym(:mio_fit_primitives_warnings), Int64,
                                        (Ptr{UInt8}, Int64), buf, n))
    (mesh=Mesh(r), num_patches=Int(rep.num_patches), num_plane=Int(rep.num_plane),
     num_sphere=Int(rep.num_sphere), num_cylinder=Int(rep.num_cylinder),
     num_cone=Int(rep.num_cone), num_torus=Int(rep.num_torus),
     num_unfitted=Int(rep.num_unfitted), num_feature_cells=Int(rep.num_feature_cells),
     warnings=isempty(text) ? String[] : String.(Base.split(text, '\n')))
end

"""
    has_occt() -> Bool

Whether the library carries the OpenCASCADE backend [`fit_cad`](@ref) needs
(a build with `-DMESHIOPLUSPLUS_WITH_OCCT=ON`).
"""
has_occt() = ccall(_sym(:mio_has_occt), Cint, ()) != 0

"""
    fit_cad(m; step_path="", tolerance=0.0, exact_edges=true,
            types=("plane", "sphere", "cylinder", "cone", "torus"),
            mesh_deflection=0.0, max_demotion_rounds=8, freeform=true,
            freeform_degree=3, freeform_max_control=64, freeform_split=true,
            fillets=true, threads=true, merge_faces=true, backend="auto",
            refine_segmentation=false, max_refine_iterations=10, refine_epsilon=1e-3,
            normal_weight=1.0, allow_merge=true, allow_split=true,
            ransac_seed=0x6d657368696f2b2b, refine_tolerance=1e-3,
            relative_tolerance=true, freeform_general=false, freeform_tolerance=0.0,
            freeform_smoothing=1e-12, max_freeform_poles=40, max_freeform_iterations=30,
            max_freeform_splits=2, sew_tolerance=0.0, extend_surfaces=false,
            snap_junctions=false, face_deviation=false, freeform_mode="legacy",
            freeform_parameterization="auto", extraordinary_vertex="g1x",
            local_refine="knots", freeform_g1_tolerance_deg=1.0,
            freeform_max_control_points=400, layout_target_patch_size=0.1,
            layout_max_valence_deviation=2, resegment_rounds=2,
            layout_snap_to_features=true, validate_by_surface_sampling=true,
            constrain_primitives=false, validation_report=false,
            validation_acceptance=false, freeform_fit_schedule="legacy",
            freeform_objective="default", kinematic_surfaces=false,
            freeform_g1_solve="sweeps")
        -> (; mesh, backend, valid, closed, tolerance, volume, skin_volume, volume_ratio,
             max_deviation, num_faces, num_plane, num_sphere, num_cylinder,
             num_cone, num_torus, num_faceted_faces, num_patches, num_demoted,
             num_exact_edges, num_approx_edges, num_solids, num_bspline,
             num_split_patches, num_fillet, num_helical, num_variable_fillet,
             num_blend, num_rounded_thread, num_vertex_blend, num_chamfer,
             num_conic_blend, num_tapered_thread, num_variable_pitch_thread,
             refined, refinement_kept, refine_iterations, num_refine_splits,
             num_refine_merges, initial_explained_area_fraction, explained_area_fraction,
             num_freeform_general, num_freeform_refinements, avg_freeform_poles,
             num_extended_edges, num_snapped_junctions, hausdorff_max, hausdorff_rms,
             brep_checked, freeform_mode, facet_area_fraction, num_freeform_adaptive,
             num_network_patches, num_network_edges, num_extraordinary_vertices,
             num_resegmented, num_escalations, max_g1_angle_deg, max_freeform_distance,
             rms_freeform_distance, num_face_records, num_constrained_primitives,
             num_constraint_groups, validation_reported, cad_to_mesh_max, cad_to_mesh_p95,
             cad_to_mesh_p99, mesh_to_cad_max, mesh_to_cad_p95, mesh_to_cad_p99,
             normal_deviation_max_deg, normal_deviation_p95_deg, step_round_trip_checked,
             step_round_trip_valid, step_round_trip_volume_ratio,
             total_freeform_control_points, num_extrusion, num_revolution,
             num_kinematic_helical, max_g1_angle_away_deg, max_g1_angle_at_ev_deg,
             num_coupled_regions)

Build a CAD solid from the fitted patches of a closed skin through the
OpenCASCADE backend: accepted fits become analytic faces, a patch no primitive
explains a rolling-ball fillet between two analytic faces, of one radius or a
radius varying along it (with `fillets`, `fit_cad:type` 7), a helical thread
face about a neighbouring cylinder's or cone's axis or the screw it fits
itself, with a straight or rounded profile (with `threads`, type 8), or a
freeform B-spline face (with `freeform`, type 6; with `freeform_split`, a patch
folding over every flattening is split into pieces first), and a patch between
two fitted faces, one of them freeform, then a rolling-ball blend (type 7).
Since v17.2.0 a patch no rolling ball or thread explains between two analytic
faces is a chamfer (type 10) or a conic blend (type 7), a thread of varying
pitch or cut on a taper is recognised, and a corner where three or more blends
meet is a vertex blend (type 9); the
rest planar facets, and a face that fails the B-rep check or
strays beyond `tolerance` is demoted to facets. `step_path` also writes the solid as STEP; the mesh is the
solid re-meshed with cell data `fit_cad:face`, `fit_cad:type` and
`fit_cad:patch`. A tolerance or deflection of 0 is automatic. Since v17.9.0,
`merge_faces` merges faces on one surface, and `backend` (`"auto"`, `"occt"`
or `"native"`) picks the engine: `"native"` needs no CAD kernel and writes the
fitted patches as analytic STEP faces, returning the skin itself; `"auto"` is
`"occt"` when the library carries it ([`has_occt`](@ref)), and `"occt"` throws
a `MeshioError` naming `-DMESHIOPLUSPLUS_WITH_OCCT=ON` without it. Since
v17.11.0, `refine_segmentation=true` first refines the patches as
[`refine_patches`](@ref) does, with `refine_tolerance` (a fraction of the
diagonal when `relative_tolerance`), `normal_weight`, `refine_epsilon`,
`max_refine_iterations`, `allow_split`, `allow_merge` and `ransac_seed`;
`refined` says whether it ran and the `refine_*`, `num_refine_*` and
`*explained_area_fraction` fields report it (zero when it did not);
`refinement_kept` says whether the solid was built from the refined patches (the
OpenCASCADE backend also builds it from the patches as given and keeps the
better). Also since
v17.11.0, `freeform_general=true` tries a patch the plane and cylinder freeform
fit leaves unexplained as a general B-spline face (`freeform_tolerance`, 0 a
quarter of the tolerance; `freeform_smoothing`; `max_freeform_poles`,
`max_freeform_iterations`, `max_freeform_splits`), reported by
`num_freeform_general`, `num_freeform_refinements` and `avg_freeform_poles`;
and, OpenCASCADE only, `sew_tolerance` (0: the tolerance) sets the sewing
tolerance, `extend_surfaces` extends freeform surfaces past their patches
(`num_extended_edges`), `snap_junctions` places the corners where analytic
faces meet by least squares (`num_snapped_junctions`), and `face_deviation`
measures each face against the skin, as the `fit_cad:face_hausdorff` and
`fit_cad:face_rms` field data and `hausdorff_max` / `hausdorff_rms` (`NaN` when
not measured). `closed` says whether every shell is closed. Since v17.12.0,
`freeform_mode` (`"legacy"`, `"adaptive"`, `"network"` or `"auto"`) picks how
the patches no primitive or feature fit explains are treated, with
`freeform_parameterization` (`"auto"`, `"projection"`, `"lscm"`,
`"harmonic"`), `extraordinary_vertex` (`"g0"`, `"g1"`, `"g1x"`, `"g1xx"`),
`local_refine` (`"knots"`, `"thb"`), `freeform_g1_tolerance_deg`,
`freeform_max_control_points`, `layout_target_patch_size`,
`layout_max_valence_deviation`, `resegment_rounds`, `layout_snap_to_features`
and `validate_by_surface_sampling`; `facet_area_fraction` is the share of the
result's area left on facet faces, `brep_checked` whether OpenCASCADE checked
the solid, and outside the legacy mode the mesh carries the
`fit_cad:face_table` field data (one row per patch tried as a surface).
`constrain_primitives` fits related primitives together on either backend
(coaxial cylinders, cones, tori and spheres share one axis; parallel and
perpendicular planes and nearly equal radii are made exact), each change kept
only where its patch still passes its fit; `num_constrained_primitives` and
`num_constraint_groups` count what it changed. `validation_report` measures
the solid against the skin both ways (`cad_to_mesh_*`, `mesh_to_cad_*`, as
maxima, 95th and 99th percentiles), the normal deviation
(`normal_deviation_max_deg`, `normal_deviation_p95_deg`) and the STEP file
read back (`step_round_trip_*`, OpenCASCADE only); `validation_acceptance`
holds every freeform face to those measures. `freeform_fit_schedule`
(`"legacy"` or `"continuation"`) and `freeform_objective` (`"default"`,
`"pdm"`, `"tdm"` or `"sdm"`) choose how the adaptive and general freeform
faces' B-spline fit runs. `kinematic_surfaces` first tries a patch no
primitive explains as a surface of linear extrusion, of revolution or a
helical surface (`num_extrusion`, `num_revolution`, `num_kinematic_helical`).
`freeform_g1_solve` (`"sweeps"` or `"coupled"`) picks how the patch network
reaches G1 (`num_coupled_regions`); `max_g1_angle_away_deg` and
`max_g1_angle_at_ev_deg` split its largest angle at the extraordinary
vertices.
See `doc/fit_cad.md`.
"""
function fit_cad(m::Mesh; step_path::AbstractString="", tolerance::Real=0.0,
                 exact_edges::Bool=true,
                 types=("plane", "sphere", "cylinder", "cone", "torus"),
                 mesh_deflection::Real=0.0, max_demotion_rounds::Integer=8,
                 freeform::Bool=true, freeform_degree::Integer=3,
                 freeform_max_control::Integer=64, freeform_split::Bool=true,
                 fillets::Bool=true, threads::Bool=true, merge_faces::Bool=true,
                 backend::AbstractString="auto", refine_segmentation::Bool=false,
                 max_refine_iterations::Integer=10, refine_epsilon::Real=1e-3,
                 normal_weight::Real=1.0, allow_merge::Bool=true, allow_split::Bool=true,
                 ransac_seed::Integer=_REFINE_DEFAULT_SEED, refine_tolerance::Real=1e-3,
                 relative_tolerance::Bool=true, freeform_general::Bool=false,
                 freeform_tolerance::Real=0.0, freeform_smoothing::Real=1e-12,
                 max_freeform_poles::Integer=40, max_freeform_iterations::Integer=30,
                 max_freeform_splits::Integer=2, sew_tolerance::Real=0.0,
                 extend_surfaces::Bool=false, snap_junctions::Bool=false,
                 face_deviation::Bool=false, freeform_mode::AbstractString="legacy",
                 freeform_parameterization::AbstractString="auto",
                 extraordinary_vertex::AbstractString="g1x",
                 local_refine::AbstractString="knots", freeform_g1_tolerance_deg::Real=1.0,
                 freeform_max_control_points::Integer=400,
                 layout_target_patch_size::Real=0.1,
                 layout_max_valence_deviation::Integer=2, resegment_rounds::Integer=2,
                 layout_snap_to_features::Bool=true,
                 validate_by_surface_sampling::Bool=true,
                 constrain_primitives::Bool=false, validation_report::Bool=false,
                 validation_acceptance::Bool=false,
                 freeform_fit_schedule::AbstractString="legacy",
                 freeform_objective::AbstractString="default",
                 kinematic_surfaces::Bool=false,
                 freeform_g1_solve::AbstractString="sweeps")
    pick(name, value, words) = begin
        k = findfirst(==(String(value)), words)
        k === nothing && throw(ArgumentError("meshio++: fit_cad: $(name) must be one of " *
                                             join(repr.(words), ", ")))
        Int32(k - 1)
    end
    mode_c = pick("freeform_mode", freeform_mode, ("legacy", "adaptive", "network", "auto"))
    param_c = pick("freeform_parameterization", freeform_parameterization,
                   ("auto", "projection", "lscm", "harmonic"))
    ev_c = pick("extraordinary_vertex", extraordinary_vertex, ("g0", "g1", "g1x", "g1xx"))
    refine_c = pick("local_refine", local_refine, ("knots", "thb"))
    schedule_c = pick("freeform_fit_schedule", freeform_fit_schedule, ("legacy", "continuation"))
    objective_c = pick("freeform_objective", freeform_objective, ("default", "pdm", "tdm", "sdm"))
    g1_solve_c = pick("freeform_g1_solve", freeform_g1_solve, ("sweeps", "coupled"))
    backend_c = backend == "auto" ? Int32(0) : backend == "occt" ? Int32(1) :
                backend == "native" ? Int32(2) :
                throw(ArgumentError("meshio++: fit_cad: backend must be \"auto\", \"occt\" or \"native\""))
    mask = UInt32(0)
    for t in types
        k = findfirst(==(String(t)), _FIT_TYPES)
        k === nothing && throw(ArgumentError("meshio++: fit_cad: unknown type $(t)"))
        mask |= UInt32(1) << k
    end
    report = Ref{_CFitCadReport}()
    report_ex = Ref{_CFitCadReportEx}()
    report_ex2 = Ref{_CFitCadReportEx2}()
    report_ex3 = Ref{_CFitCadReportEx3}()
    report_ex4 = Ref{_CFitCadReportEx4}()
    opts = _CFitCadOpts(Float64(tolerance), Float64(mesh_deflection),
                        Int64(max_demotion_rounds), Int32(exact_edges), mask,
                        Int32(freeform), Int32(freeform_degree), Int64(freeform_max_control),
                        Int32(freeform_split), Int32(fillets), Int32(threads),
                        Int32(merge_faces), backend_c, Int32(0), ntuple(_ -> Int64(0), 1))
    # The refinement's other options keep their C defaults; a NULL prefix is "patch".
    refine = _CRefinePatchesOpts(Float64(refine_tolerance), Float64(normal_weight),
                                 Float64(refine_epsilon), 30.0, Int64(max_refine_iterations),
                                 Int64(3), UInt64(ransac_seed), _FIT_ALL_TYPES,
                                 Int32(relative_tolerance), Int32(allow_split),
                                 Int32(allow_merge), Cstring(C_NULL), ntuple(_ -> Int64(0), 8))
    opts_ex = _CFitCadOptsEx(Int32(refine_segmentation), Int32(freeform_general), refine,
                             Float64(freeform_tolerance), Float64(freeform_smoothing),
                             Int64(max_freeform_poles), Int64(max_freeform_iterations),
                             Int64(max_freeform_splits), Float64(sew_tolerance),
                             Int32(extend_surfaces), Int32(snap_junctions),
                             Int32(face_deviation), Int32(0), ntuple(_ -> Int64(0), 8))
    opts_ex2 = _CFitCadOptsEx2(mode_c, param_c, ev_c, refine_c,
                               Float64(freeform_g1_tolerance_deg),
                               Int64(freeform_max_control_points),
                               Float64(layout_target_patch_size),
                               Int64(layout_max_valence_deviation), Int64(resegment_rounds),
                               Int32(layout_snap_to_features),
                               Int32(validate_by_surface_sampling),
                               Int32(constrain_primitives), Int32(validation_report),
                               Int32(validation_acceptance), schedule_c, objective_c,
                               Int32(kinematic_surfaces), g1_solve_c, Int32(0),
                               ntuple(_ -> Int64(0), 4))
    ptr = ccall(_sym(:mio_fit_cad4), Ptr{Cvoid},
                (Ptr{Cvoid}, Cstring, Ref{_CFitCadOpts}, Ref{_CFitCadOptsEx},
                 Ref{_CFitCadOptsEx2}, Ptr{_CFitCadReport}, Ptr{_CFitCadReportEx},
                 Ptr{_CFitCadReportEx2}, Ptr{_CFitCadReportEx3}, Ptr{_CFitCadReportEx4}),
                _handle(m), String(step_path), Ref(opts), Ref(opts_ex), Ref(opts_ex2), report,
                report_ex, report_ex2, report_ex3, report_ex4)
    r = _check_ptr(ptr)
    rep = report[]
    ex = report_ex[]
    ex2 = report_ex2[]
    ex3 = report_ex3[]
    ex4 = report_ex4[]
    (mesh=Mesh(r), backend=ex.backend == 2 ? "native" : "occt", valid=rep.valid != 0,
     closed=ex2.closed != 0, tolerance=rep.tolerance, volume=rep.volume,
     skin_volume=rep.skin_volume, volume_ratio=rep.volume_ratio,
     max_deviation=rep.max_deviation, num_faces=Int(rep.num_faces),
     num_plane=Int(rep.num_plane), num_sphere=Int(rep.num_sphere),
     num_cylinder=Int(rep.num_cylinder), num_cone=Int(rep.num_cone),
     num_torus=Int(rep.num_torus), num_faceted_faces=Int(rep.num_faceted_faces),
     num_patches=Int(rep.num_patches), num_demoted=Int(rep.num_demoted),
     num_exact_edges=Int(rep.num_exact_edges), num_approx_edges=Int(rep.num_approx_edges),
     num_solids=Int(rep.num_solids), num_bspline=Int(rep.num_bspline),
     num_split_patches=Int(rep.num_split_patches), num_fillet=Int(rep.num_fillet),
     num_helical=Int(rep.num_helical), num_variable_fillet=Int(rep.num_variable_fillet),
     num_blend=Int(rep.num_blend), num_rounded_thread=Int(rep.num_rounded_thread),
     num_vertex_blend=Int(ex.num_vertex_blend), num_chamfer=Int(ex.num_chamfer),
     num_conic_blend=Int(ex.num_conic_blend), num_tapered_thread=Int(ex.num_tapered_thread),
     num_variable_pitch_thread=Int(ex.num_variable_pitch_thread),
     refined=ex2.refined != 0, refinement_kept=ex2.refinement_kept != 0,
     refine_iterations=Int(ex2.refine_iterations),
     num_refine_splits=Int(ex2.num_refine_splits),
     num_refine_merges=Int(ex2.num_refine_merges),
     initial_explained_area_fraction=ex2.initial_explained_area_fraction,
     explained_area_fraction=ex2.explained_area_fraction,
     num_freeform_general=Int(ex2.num_freeform_general),
     num_freeform_refinements=Int(ex2.num_freeform_refinements),
     avg_freeform_poles=ex2.avg_freeform_poles,
     num_extended_edges=Int(ex2.num_extended_edges),
     num_snapped_junctions=Int(ex2.num_snapped_junctions),
     hausdorff_max=ex2.hausdorff_max, hausdorff_rms=ex2.hausdorff_rms,
     brep_checked=ex3.brep_checked != 0,
     freeform_mode=("legacy", "adaptive", "network", "auto")[ex3.freeform_mode + 1],
     facet_area_fraction=ex3.facet_area_fraction,
     num_freeform_adaptive=Int(ex3.num_freeform_adaptive),
     num_network_patches=Int(ex3.num_network_patches),
     num_network_edges=Int(ex3.num_network_edges),
     num_extraordinary_vertices=Int(ex3.num_extraordinary_vertices),
     num_resegmented=Int(ex3.num_resegmented), num_escalations=Int(ex3.num_escalations),
     max_g1_angle_deg=ex3.max_g1_angle_deg, max_freeform_distance=ex3.max_freeform_distance,
     rms_freeform_distance=ex3.rms_freeform_distance,
     num_face_records=Int(ex3.num_face_records),
     num_constrained_primitives=Int(ex3.num_constrained_primitives),
     num_constraint_groups=Int(ex3.num_constraint_groups),
     validation_reported=ex4.validation_reported != 0,
     cad_to_mesh_max=ex4.cad_to_mesh_max, cad_to_mesh_p95=ex4.cad_to_mesh_p95,
     cad_to_mesh_p99=ex4.cad_to_mesh_p99, mesh_to_cad_max=ex4.mesh_to_cad_max,
     mesh_to_cad_p95=ex4.mesh_to_cad_p95, mesh_to_cad_p99=ex4.mesh_to_cad_p99,
     normal_deviation_max_deg=ex4.normal_deviation_max_deg,
     normal_deviation_p95_deg=ex4.normal_deviation_p95_deg,
     step_round_trip_checked=ex4.step_round_trip_checked != 0,
     step_round_trip_valid=ex4.step_round_trip_valid != 0,
     step_round_trip_volume_ratio=ex4.step_round_trip_volume_ratio,
     total_freeform_control_points=Int(ex4.total_freeform_control_points),
     num_extrusion=Int(ex4.num_extrusion), num_revolution=Int(ex4.num_revolution),
     num_kinematic_helical=Int(ex4.num_kinematic_helical),
     max_g1_angle_away_deg=ex4.max_g1_angle_away_deg,
     max_g1_angle_at_ev_deg=ex4.max_g1_angle_at_ev_deg,
     num_coupled_regions=Int(ex4.num_coupled_regions))
end

"""
    refine_patches(m; tolerance=1e-3, relative=true, normal_weight=1.0, epsilon=1e-3,
                   feature_angle=30.0, max_iterations=10, min_patch_cells=3,
                   ransac_seed=0x6d657368696f2b2b, types=nothing, allow_split=true,
                   allow_merge=true, prefix="patch")
        -> (; mesh, num_patches, num_iterations, num_splits, num_merges, num_absorbed,
             num_primitive, num_freeform, num_facets, tolerance, explained_area_fraction,
             initial_explained_area_fraction, hausdorff_max, hausdorff_rms,
             initial_energy, energy, stop)

Refine the `segment:patch` labels of a surface (from [`segment`](@ref)) by
variational shape approximation (v17.11.0): patch boundaries move, a freeform
patch is split where seeded RANSAC finds a primitive inside it, and two
neighbouring patches one primitive explains are merged, until each patch is
one plane, sphere, cylinder, cone or torus, or a freeform proxy. `tolerance` is
a fraction of the bounding-box diagonal when `relative`; `types` (a collection
of `"plane"`, `"sphere"`, `"cylinder"`, `"cone"`, `"torus"`; `nothing` for
all) limits the proxies; `ransac_seed` is a `UInt64`, exact over its whole
range. The mesh comes back with `segment:patch`, `segment:adjacency` and the
`<prefix>:<k>` regions rewritten, the `fit:*` data of
[`fit_primitives`](@ref) on the refined patches and per-patch `refine:*`
field data; `stop` is `:unchanged`, `:small_gain`, `:iteration_cap` or
`:energy_rose`. Throws a `MeshioError` without `segment:patch` or on an
option out of range. See `doc/refine_patches.md`.
"""
function refine_patches(m::Mesh; tolerance::Real=1e-3, relative::Bool=true,
                        normal_weight::Real=1.0, epsilon::Real=1e-3,
                        feature_angle::Real=30.0, max_iterations::Integer=10,
                        min_patch_cells::Integer=3,
                        ransac_seed::Integer=_REFINE_DEFAULT_SEED, types=nothing,
                        allow_split::Bool=true, allow_merge::Bool=true,
                        prefix::AbstractString="patch")
    mask = types === nothing ? _FIT_ALL_TYPES : _fit_type_mask(types, "refine_patches")
    report = Ref{_CRefinePatchesReport}()
    prefix_c = Vector{UInt8}(codeunits(String(prefix) * "\0"))
    ptr = GC.@preserve prefix_c begin
        opts = _CRefinePatchesOpts(Float64(tolerance), Float64(normal_weight),
                                   Float64(epsilon), Float64(feature_angle),
                                   Int64(max_iterations), Int64(min_patch_cells),
                                   UInt64(ransac_seed), mask, Int32(relative),
                                   Int32(allow_split), Int32(allow_merge),
                                   Cstring(pointer(prefix_c)), ntuple(_ -> Int64(0), 8))
        ccall(_sym(:mio_refine_patches), Ptr{Cvoid},
              (Ptr{Cvoid}, Ref{_CRefinePatchesOpts}, Ptr{_CRefinePatchesReport}),
              _handle(m), Ref(opts), report)
    end
    r = _check_ptr(ptr)
    rep = report[]
    stop = 0 <= rep.stop < length(_REFINE_STOPS) ? _REFINE_STOPS[rep.stop + 1] : :unknown
    (mesh=Mesh(r), num_patches=Int(rep.num_patches), num_iterations=Int(rep.num_iterations),
     num_splits=Int(rep.num_splits), num_merges=Int(rep.num_merges),
     num_absorbed=Int(rep.num_absorbed), num_primitive=Int(rep.num_primitive),
     num_freeform=Int(rep.num_freeform), num_facets=Int(rep.num_facets),
     tolerance=rep.tolerance, explained_area_fraction=rep.explained_area_fraction,
     initial_explained_area_fraction=rep.initial_explained_area_fraction,
     hausdorff_max=rep.hausdorff_max, hausdorff_rms=rep.hausdorff_rms,
     initial_energy=rep.initial_energy, energy=rep.energy, stop=stop)
end

"""
    hausdorff_distance(a, b; face_samples=0, region_a="", region_b="")
        -> (; distance, a_to_b, b_to_a, mean_a_to_b, rms_a_to_b, mean_b_to_a,
             rms_b_to_a, num_samples_a, num_samples_b, worst_point_a, worst_point_b)

The (sampled) Hausdorff distance between the surfaces of two meshes; a volume
mesh contributes its skin. `face_samples = s > 0` also samples the centroids of
the `s*s` sub-triangles of every triangle -- vertex sampling alone is a lower
bound. See `doc/hausdorff.md`.
"""
function hausdorff_distance(a::Mesh, b::Mesh; face_samples::Integer=0,
                            region_a::AbstractString="", region_b::AbstractString="")
    report = Ref{_CHausdorffReport}()
    ra = Vector{UInt8}(codeunits(String(region_a) * "\0"))
    rb = Vector{UInt8}(codeunits(String(region_b) * "\0"))
    GC.@preserve ra rb begin
        opts = _CHausdorffOpts(Cstring(pointer(ra)), Cstring(pointer(rb)), Int64(face_samples),
                               0.0, ntuple(_ -> Int64(0), 6))
        _check(ccall(_sym(:mio_hausdorff_distance), Cint,
                     (Ptr{Cvoid}, Ptr{Cvoid}, Ref{_CHausdorffOpts}, Ptr{_CHausdorffReport}),
                     _handle(a), _handle(b), Ref(opts), report))
    end
    rep = report[]
    (distance=rep.distance, a_to_b=rep.a_to_b, b_to_a=rep.b_to_a,
     mean_a_to_b=rep.mean_a_to_b, rms_a_to_b=rep.rms_a_to_b, mean_b_to_a=rep.mean_b_to_a,
     rms_b_to_a=rep.rms_b_to_a, num_samples_a=Int(rep.num_samples_a),
     num_samples_b=Int(rep.num_samples_b), worst_point_a=collect(rep.worst_point_a),
     worst_point_b=collect(rep.worst_point_b))
end

"""
    check_quality(m; require="", max_inverted=0, max_degenerate=0)
        -> (; passed, num_checks, num_failed, num_cells, num_inverted, num_degenerate, summary)

Quality gate: score the cells and test the `require` thresholds (text such as
`"scaled_jacobian >= 0.2; aspect_ratio <= 5 @ 1%"`) plus the inverted and
degenerate counts (negative disables). `summary` is the text both CLIs print.
See `doc/mesh_quality.md`.
"""
function check_quality(m::Mesh; require::AbstractString="", max_inverted::Integer=0,
                       max_degenerate::Integer=0)
    report = Ref{_CQualityGateReport}()
    buf = zeros(UInt8, 8192)
    _check(ccall(_sym(:mio_check_quality), Cint,
                 (Ptr{Cvoid}, Cstring, Int64, Int64, Ptr{_CQualityGateReport}, Ptr{UInt8}, Int64),
                 _handle(m), String(require), Int64(max_inverted), Int64(max_degenerate),
                 report, buf, Int64(length(buf))))
    r = report[]
    n = min(Int(r.summary_length), length(buf) - 1)
    (passed=r.passed != 0, num_checks=Int(r.num_checks), num_failed=Int(r.num_failed),
     num_cells=Int(r.num_cells), num_inverted=Int(r.num_inverted),
     num_degenerate=Int(r.num_degenerate), summary=String(buf[1:n]))
end

_quality_tuple(q::_CSurfaceQuality) =
    (boundary_edges=Int(q.boundary_edges), non_manifold_edges=Int(q.non_manifold_edges),
     inconsistent_pairs=Int(q.inconsistent_pairs),
     degenerate_triangles=Int(q.degenerate_triangles), watertight=q.watertight != 0)

"""
    repair(m; fix_orientation=true, orient_outward=true, fill_holes=true,
           split_non_manifold=true, max_hole_edges=10, weld_tolerance=0.0, fill_style="fan",
           record_provenance=false)
        -> (; mesh, quality_before, quality_after, num_flipped, num_components,
             largest_component, num_oriented_outward, num_unorientable,
             num_vertices_split, num_holes_detected, num_holes_filled,
             num_holes_skipped, num_faces_added, num_points_added, points_welded)

Repair a surface's orientation, holes and pinched vertices: weld (opt-in) ->
triangulate (blocks 1:1) -> split bowties -> orient by the topological
half-edge rule per connected component -> fan-fill boundary loops of at most
`max_hole_edges` edges (`<= 0` means no limit) -> orient closed components
outward. Lower-dimensional blocks ride along; fill triangles land in one
trailing `triangle` block. The C++ index maps are not exposed on the flat ABI.
See `doc/repair.md`.
"""
function repair(m::Mesh; fix_orientation::Bool=true, orient_outward::Bool=true,
                fill_holes::Bool=true, split_non_manifold::Bool=true,
                max_hole_edges::Integer=10, weld_tolerance::Real=0.0,
                record_provenance::Bool=false, fill_style::AbstractString="fan")
    style = fill_style == "fan" ? Int32(0) : fill_style == "flat" ? Int32(1) :
            fill_style == "curved" ? Int32(2) :
            throw(ArgumentError("meshio++: repair: fill_style must be fan, flat or curved"))
    report = Ref{_CRepairReport}()
    opts = _CRepairOpts(fix_orientation ? Int32(1) : Int32(0),
                        orient_outward ? Int32(1) : Int32(0),
                        fill_holes ? Int32(1) : Int32(0),
                        split_non_manifold ? Int32(1) : Int32(0),
                        record_provenance ? Int32(1) : Int32(0), style,
                        Int64(max_hole_edges), Float64(weld_tolerance),
                        (Int64(0), Int64(0), Int64(0), Int64(0), Int64(0)))
    ptr = ccall(_sym(:mio_repair), Ptr{Cvoid},
                (Ptr{Cvoid}, Ref{_CRepairOpts}, Ptr{_CRepairReport}),
                _handle(m), Ref(opts), report)
    r = _check_ptr(ptr)
    rep = report[]
    (mesh=Mesh(r), quality_before=_quality_tuple(rep.quality_before),
     quality_after=_quality_tuple(rep.quality_after),
     num_flipped=Int(rep.num_flipped), num_components=Int(rep.num_components),
     largest_component=Int(rep.largest_component),
     num_oriented_outward=Int(rep.num_oriented_outward),
     num_unorientable=Int(rep.num_unorientable),
     num_vertices_split=Int(rep.num_vertices_split),
     num_holes_detected=Int(rep.num_holes_detected),
     num_holes_filled=Int(rep.num_holes_filled),
     num_holes_skipped=Int(rep.num_holes_skipped),
     num_faces_added=Int(rep.num_faces_added),
     num_points_added=Int(rep.num_points_added), points_welded=Int(rep.points_welded))
end

const _SHRINKWRAP_WEIGHTS = Dict(:angle => Int32(0), :area => Int32(1))

"""
    shrinkwrap(m, target; offset=0.0, max_distance=0.0, weights="",
               target_region="", normal_weight=:angle, record_distance=false,
               record_closest_cell=false)
        -> (; mesh, quality, num_projected, num_missed, num_skipped, max_displacement)

Project every (selected) point of `m` onto the surface of `target`:
`x' = x + w (p + offset n - x)`, one projection, no iteration. Every point of
the source moves whatever cells it carries; only the target must be a surface.
`weights` names a point-data array on the source (a float array blends, an
integer/bool one selects); the offset goes along the hit FEATURE's
pseudonormal (the bisector at a crease). `quality` is the target's. See
`doc/shrinkwrap.md`.
"""
function shrinkwrap(m::Mesh, target::Mesh; offset::Real=0.0, max_distance::Real=0.0,
                    weights::AbstractString="", target_region::AbstractString="",
                    normal_weight::Symbol=:angle, record_distance::Bool=false,
                    record_closest_cell::Bool=false)
    haskey(_SHRINKWRAP_WEIGHTS, normal_weight) ||
        throw(ArgumentError("meshio++: shrinkwrap: unknown normal weight '$(normal_weight)' " *
                            "(expected :angle or :area)"))
    report = Ref{_CShrinkwrapReport}()
    weights_c = Vector{UInt8}(codeunits(String(weights) * "\0"))
    region_c = Vector{UInt8}(codeunits(String(target_region) * "\0"))
    ptr = GC.@preserve weights_c region_c begin
        opts = _CShrinkwrapOpts(Cstring(pointer(weights_c)), Cstring(pointer(region_c)),
                                Float64(offset), Float64(max_distance), 0.0,
                                _SHRINKWRAP_WEIGHTS[normal_weight],
                                record_distance ? Int32(1) : Int32(0),
                                record_closest_cell ? Int32(1) : Int32(0), Int32(0),
                                (Int64(0), Int64(0), Int64(0), Int64(0), Int64(0)))
        ccall(_sym(:mio_shrinkwrap), Ptr{Cvoid},
              (Ptr{Cvoid}, Ptr{Cvoid}, Ref{_CShrinkwrapOpts}, Ptr{_CShrinkwrapReport}),
              _handle(m), _handle(target), Ref(opts), report)
    end
    r = _check_ptr(ptr)
    rep = report[]
    (mesh=Mesh(r), quality=_quality_tuple(rep.quality),
     num_projected=Int(rep.num_projected), num_missed=Int(rep.num_missed),
     num_skipped=Int(rep.num_skipped), max_displacement=Float64(rep.max_displacement))
end

"""
    sobolev_deform(m, array, length_scale; fixed_points_array="",
                   fix_boundary=false, record_filtered=false,
                   max_iterations=128, tolerance=1e-10)
        -> (; mesh, num_iterations, residual, converged, num_fixed, num_isolated,
             max_displacement)

Sobolev (Helmholtz-filtered) deformation: solve `(M + l^2 K) u = M d` over the
mesh's own P1 operators (`K` from the simplex edge Gram matrix, `M` the uniform
mean lumped mass, matrix-free Jacobi-PCG) and move the points by `u`. `array`
names the `(n, dim)` point-data displacement; `fixed_points_array` an
integer/bool point-data array whose nonzero entries pin their point (the C++
mask form is a flat-ABI gap). Every top-dimensional block must be a linear
simplex; lower-dimensional blocks ride along. Non-convergence is reported
through `converged`, never thrown. See `doc/sobolev_deform.md`.
"""
function sobolev_deform(m::Mesh, array::AbstractString, length_scale::Real;
                        fixed_points_array::AbstractString="", fix_boundary::Bool=false,
                        record_filtered::Bool=false, max_iterations::Integer=128,
                        tolerance::Real=1e-10)
    report = Ref{_CSobolevReport}()
    array_c = Vector{UInt8}(codeunits(String(array) * "\0"))
    fixed_c = Vector{UInt8}(codeunits(String(fixed_points_array) * "\0"))
    ptr = GC.@preserve array_c fixed_c begin
        opts = _CSobolevOpts(Cstring(pointer(array_c)), Cstring(pointer(fixed_c)),
                             Float64(length_scale), Float64(tolerance),
                             Int32(max_iterations), fix_boundary ? Int32(1) : Int32(0),
                             record_filtered ? Int32(1) : Int32(0), Int32(0),
                             (Int64(0), Int64(0), Int64(0), Int64(0), Int64(0)))
        ccall(_sym(:mio_sobolev_deform), Ptr{Cvoid},
              (Ptr{Cvoid}, Ref{_CSobolevOpts}, Ptr{_CSobolevReport}),
              _handle(m), Ref(opts), report)
    end
    r = _check_ptr(ptr)
    rep = report[]
    (mesh=Mesh(r), num_iterations=Int(rep.num_iterations), residual=Float64(rep.residual),
     converged=rep.converged != 0, num_fixed=Int(rep.num_fixed),
     num_isolated=Int(rep.num_isolated), max_displacement=Float64(rep.max_displacement))
end

"""
    prepare_metric(m, arrays; output="metric", target_vertices=0.0, norm=2.0, hmin=0.0,
                   hmax=0.0, gradation=0.0, max_sweeps=100)
        -> (; mesh, complexity_before, complexity_after, num_graded, sweeps)

The metric toolbox: intersect the metric fields `arrays`, normalise to
`target_vertices`, clamp to `[hmin, hmax]` and grade, into point data `output`.
See `doc/metric.md`.
"""
function prepare_metric(m::Mesh, arrays::AbstractVector{<:AbstractString}; output::AbstractString="metric",
                        target_vertices::Real=0.0, norm::Real=2.0, hmin::Real=0.0, hmax::Real=0.0,
                        gradation::Real=0.0, max_sweeps::Integer=100)
    names = [Vector{UInt8}(codeunits(String(a) * "\0")) for a in arrays]
    ptrs = [Cstring(pointer(n)) for n in names]
    rep = Ref{_CMetricReport}()
    ptr = GC.@preserve names ptrs output begin
        opts = _CMetricOpts(isempty(ptrs) ? Ptr{Cstring}(C_NULL) : pointer(ptrs), Int64(length(ptrs)),
                            Base.unsafe_convert(Cstring, output), Float64(target_vertices),
                            Float64(norm), Float64(hmin), Float64(hmax), Float64(gradation),
                            Int64(max_sweeps), ntuple(_ -> Int64(0), 4))
        ccall(_sym(:mio_prepare_metric), Ptr{Cvoid},
              (Ptr{Cvoid}, Ref{_CMetricOpts}, Ref{_CMetricReport}), _handle(m), Ref(opts), rep)
    end
    r = _check_ptr(ptr)
    h = rep[]
    (mesh=Mesh(r), complexity_before=h.complexity_before, complexity_after=h.complexity_after,
     num_graded=Int(h.num_graded), sweeps=Int(h.sweeps))
end

"""
    untangle(m; fix_boundary=true, rings=2, smooth=false, max_sweeps=100)
        -> (; mesh, num_inverted_before, num_inverted_after, min_quality_before,
              min_quality_after, sweeps, moves)

Untangle tetrahedra and hexahedra: move nodes until every corner has a positive
signed volume (Escobar et al.'s functional). See `doc/untangle.md`.
"""
function untangle(m::Mesh; fix_boundary::Bool=true, rings::Integer=2, smooth::Bool=false,
                  max_sweeps::Integer=100)
    opts = _CUntangleOpts(Ptr{UInt8}(C_NULL), Int64(0), Int64(rings), Int64(max_sweeps),
                          fix_boundary ? Int32(1) : Int32(0), smooth ? Int32(1) : Int32(0),
                          ntuple(_ -> Int64(0), 4))
    rep = Ref{_CUntangleReport}()
    ptr = ccall(_sym(:mio_untangle), Ptr{Cvoid},
                (Ptr{Cvoid}, Ref{_CUntangleOpts}, Ref{_CUntangleReport}), _handle(m), Ref(opts), rep)
    r = _check_ptr(ptr)
    h = rep[]
    (mesh=Mesh(r), num_inverted_before=Int(h.num_inverted_before),
     num_inverted_after=Int(h.num_inverted_after), min_quality_before=h.min_quality_before,
     min_quality_after=h.min_quality_after, sweeps=Int(h.sweeps), moves=Int(h.moves))
end

"""
    morph(m, array; fix_boundary=true, constrained_array="", stiffening=1.0, steps=1,
          allow_inverted=false, tolerance=1e-10, max_iterations=0)
        -> (; mesh, num_constrained, num_free, num_inverted, iterations, max_displacement, converged)

Boundary-driven morphing: the constrained points (the boundary, and the
nonzero entries of point data `constrained_array`) move by point data `array`,
the rest by its harmonic extension with Jacobian-based stiffening; inverted
cells raise unless `allow_inverted`. See `doc/morph.md`.
"""
function morph(m::Mesh, array::AbstractString; fix_boundary::Bool=true,
               constrained_array::AbstractString="", stiffening::Real=1.0, steps::Integer=1,
               allow_inverted::Bool=false, tolerance::Real=1e-10, max_iterations::Integer=0)
    rep = Ref{_CMorphReport}()
    ptr = GC.@preserve array constrained_array begin
        opts = _CMorphOpts(Base.unsafe_convert(Cstring, array),
                           isempty(constrained_array) ? Cstring(C_NULL) :
                           Base.unsafe_convert(Cstring, constrained_array),
                           Ptr{UInt8}(C_NULL), Int64(0), Float64(stiffening), Int64(steps),
                           Float64(tolerance), Int64(max_iterations),
                           fix_boundary ? Int32(1) : Int32(0), allow_inverted ? Int32(1) : Int32(0),
                           ntuple(_ -> Int64(0), 4))
        ccall(_sym(:mio_morph), Ptr{Cvoid}, (Ptr{Cvoid}, Ref{_CMorphOpts}, Ref{_CMorphReport}),
              _handle(m), Ref(opts), rep)
    end
    r = _check_ptr(ptr)
    h = rep[]
    (mesh=Mesh(r), num_constrained=Int(h.num_constrained), num_free=Int(h.num_free),
     num_inverted=Int(h.num_inverted), iterations=Int(h.iterations),
     max_displacement=h.max_displacement, converged=h.converged != 0)
end

"""
    laplacian_eigenmodes(m; num_modes=10, tolerance=1e-8, max_iterations=0)
        -> (; mesh, iterations, max_residual, converged, num_degenerate)

The lowest eigenpairs of a triangle surface's Laplace-Beltrami operator: the
mesh carries point data `laplacian:modes` (n x num_modes, M-orthonormal) and
field data `laplacian:eigenvalues` (ascending). See `doc/eigenmodes.md`.
"""
function laplacian_eigenmodes(m::Mesh; num_modes::Integer=10, tolerance::Real=1e-8,
                              max_iterations::Integer=0)
    opts = _CEigenmodeOpts(Int64(num_modes), Float64(tolerance), Int64(max_iterations),
                           ntuple(_ -> Int64(0), 4))
    rep = Ref{_CEigenmodeReport}()
    ptr = ccall(_sym(:mio_laplacian_eigenmodes), Ptr{Cvoid},
                (Ptr{Cvoid}, Ref{_CEigenmodeOpts}, Ref{_CEigenmodeReport}),
                _handle(m), Ref(opts), rep)
    r = _check_ptr(ptr)
    h = rep[]
    (mesh=Mesh(r), iterations=Int(h.iterations), max_residual=h.max_residual,
     converged=h.converged != 0, num_degenerate=Int(h.num_degenerate))
end

"""
    register_rigid(source, target; method="point_to_point", max_iterations=50,
                   tolerance=1e-10, max_distance=0.0, initial=nothing,
                   symmetric_check=false, symmetric_tolerance=0.05, source_region="",
                   target_region="", adaptive_trim=0.0)
        -> (; matrix, initial_rms, rms, iterations, num_pairs, converged, num_rejected,
              candidate_rms, chosen_candidate, warnings)

Rigid registration by iterative closest points: `matrix` (4x4, `matrix[i, j]`
row `i`) lays `source`'s points onto `target` (its surface when it has surface
cells, else its points); apply it with `transform(source; matrix=...)`.
`initial` is a starting 4x4 matrix, or `"auto"`: the best of the 24
principal-frame candidates (roadmap 7.3.1; `chosen_candidate` is 0-based).
`source_region`/`target_region` restrict the pairs to cell regions and
`adaptive_trim` drops pairs beyond that multiple of the median pair distance
(7.3.3). See `doc/registration.md`.
"""
function register_rigid(source::Mesh, target::Mesh; method::AbstractString="point_to_point",
                        max_iterations::Integer=50, tolerance::Real=1e-10,
                        max_distance::Real=0.0, initial=nothing, symmetric_check::Bool=false,
                        symmetric_tolerance::Real=0.05, source_region::AbstractString="",
                        target_region::AbstractString="", adaptive_trim::Real=0.0)
    meth = method == "point_to_point" ? Int32(0) : method == "point_to_plane" ? Int32(1) :
           throw(ArgumentError("meshio++: register_rigid: method must be \"point_to_point\" or \"point_to_plane\""))
    mode = Int32(0)
    if initial isa AbstractString
        initial == "auto" ||
            throw(ArgumentError("meshio++: register_rigid: initial must be a 4x4 matrix, nothing or \"auto\""))
        mode = Int32(1)
        initial = nothing
    end
    start = initial === nothing ? Matrix{Float64}(I4) : Matrix{Float64}(initial)
    size(start) == (4, 4) || throw(ArgumentError("meshio++: register_rigid: initial must be 4x4"))
    flat = ntuple(k -> start[(k - 1) ÷ 4 + 1, (k - 1) % 4 + 1], 16)
    base = _CRegistrationOpts(meth, Int32(0), Int64(max_iterations), Float64(tolerance),
                              Float64(max_distance), flat, ntuple(_ -> Int64(0), 4))
    rep = Ref{_CRegistrationReport}()
    sr = Vector{UInt8}(codeunits(String(source_region) * "\0"))
    tr = Vector{UInt8}(codeunits(String(target_region) * "\0"))
    st = GC.@preserve sr tr begin
        opts = _CRegistrationExOpts(base, mode, symmetric_check ? Int32(1) : Int32(0),
                                    Cdouble(symmetric_tolerance), Cdouble(adaptive_trim),
                                    Cstring(pointer(sr)), Cstring(pointer(tr)),
                                    ntuple(_ -> Int64(0), 4))
        ccall(_sym(:mio_register_rigid_ex), Cint,
              (Ptr{Cvoid}, Ptr{Cvoid}, Ref{_CRegistrationExOpts}, Ref{_CRegistrationReport}),
              _handle(source), _handle(target), Ref(opts), rep)
    end
    _check(st)
    r = rep[]
    matrix = [r.matrix[4 * (i - 1) + j] for i in 1:4, j in 1:4]
    cands = zeros(Cdouble, Int(r.num_candidates))
    ccall(_sym(:mio_register_rigid_candidates), Int64, (Ptr{Cdouble}, Int64), cands,
          length(cands))
    text = _getstring((buf, n) -> ccall(_sym(:mio_register_rigid_warnings), Int64,
                                        (Ptr{UInt8}, Int64), buf, n))
    (matrix=matrix, initial_rms=r.initial_rms, rms=r.rms, iterations=Int(r.iterations),
     num_pairs=Int(r.num_pairs), converged=r.converged != 0,
     num_rejected=Int(r.num_rejected), candidate_rms=cands,
     chosen_candidate=Int(r.chosen_candidate),
     warnings=isempty(text) ? String[] : String.(Base.split(text, '\n')))
end

"""
    register_landmarks(source_points, target_points; weights=nothing)
        -> (; matrix, residuals, rms)

The rigid motion laying picked points (rows of `n x 3` matrices, `n >= 3`, not
collinear) onto their targets: the weighted Horn solve (roadmap 7.3.2).
"""
function register_landmarks(source_points::AbstractMatrix{<:Real},
                            target_points::AbstractMatrix{<:Real}; weights=nothing)
    n = size(source_points, 1)
    (size(source_points, 2) == 3 && size(target_points) == (n, 3)) ||
        throw(ArgumentError("meshio++: register_landmarks: points must be n x 3 matrices"))
    src = Cdouble.(vec(permutedims(source_points)))
    dst = Cdouble.(vec(permutedims(target_points)))
    w = weights === nothing ? C_NULL : Cdouble.(collect(weights))
    m = zeros(Cdouble, 16)
    res = zeros(Cdouble, max(n, 1))
    rms = Ref{Cdouble}(0.0)
    st = ccall(_sym(:mio_register_landmarks), Cint,
               (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Int64, Ptr{Cdouble}, Ptr{Cdouble},
                Ref{Cdouble}), src, dst, w, n, m, res, rms)
    _check(st)
    (matrix=[m[4 * (i - 1) + j] for i in 1:4, j in 1:4], residuals=res[1:n], rms=rms[])
end

"""
    find_symmetry_plane(m; seed=nothing, tolerance=0.0, min_share=0.9, max_samples=4000,
                        max_iterations=50) -> (; origin, normal, rms, share, tolerance,
                                                chosen_candidate)

The plane `m` is most nearly symmetric about (roadmap 7.3.6): each principal
plane (and `seed`, an `(origin, normal)` pair) is mirrored, registered back by
ICP and read off the reflection-motion product. Refuses a plane below
`min_share`. See `doc/registration.md`.
"""
function find_symmetry_plane(m::Mesh; seed=nothing, tolerance::Real=0.0,
                             min_share::Real=0.9, max_samples::Integer=4000,
                             max_iterations::Integer=50)
    so, sn = seed === nothing ? ((0.0, 0.0, 0.0), (0.0, 0.0, 0.0)) :
             (Tuple(Float64.(seed[1])), Tuple(Float64.(seed[2])))
    opts = _CSymmetryPlaneOpts(so, sn, Cdouble(tolerance), Cdouble(min_share),
                               Int64(max_samples), Int64(max_iterations), ntuple(_ -> Int64(0), 4))
    rep = Ref{_CSymmetryPlaneReport}()
    st = ccall(_sym(:mio_find_symmetry_plane), Cint,
               (Ptr{Cvoid}, Ref{_CSymmetryPlaneOpts}, Ref{_CSymmetryPlaneReport}),
               _handle(m), Ref(opts), rep)
    _check(st)
    r = rep[]
    (origin=collect(r.origin), normal=collect(r.normal), rms=r.rms, share=r.share,
     tolerance=r.tolerance, chosen_candidate=Int(r.chosen_candidate))
end

"""
    auto_align(m; angle_tolerance=10.0, subdivisions=3, region="", rotate_vector_data=false)
        -> (; mesh, matrix, mode, primary, secondary, share, face_share, axis_share,
              opening_area)

Pose a scan from itself alone (roadmap 7.3.5): on its largest flat face
(`mode == "face"`) or a cylinder on the end of its axis (`"axis"`), from an
area-weighted histogram of its facet normals. See `doc/registration.md`.
"""
function auto_align(m::Mesh; angle_tolerance::Real=10.0, subdivisions::Integer=3,
                    region::AbstractString="", rotate_vector_data::Bool=false)
    rg = Vector{UInt8}(codeunits(String(region) * "\0"))
    rep = Ref{_CAutoAlignReport}()
    ptr = GC.@preserve rg begin
        opts = _CAutoAlignOpts(Cdouble(angle_tolerance), Int64(subdivisions), Cstring(pointer(rg)),
                               rotate_vector_data ? Int32(1) : Int32(0), Int32(0),
                               ntuple(_ -> Int64(0), 4))
        ccall(_sym(:mio_auto_align), Ptr{Cvoid},
              (Ptr{Cvoid}, Ref{_CAutoAlignOpts}, Ref{_CAutoAlignReport}), _handle(m), Ref(opts), rep)
    end
    out = Mesh(_check_ptr(ptr))
    r = rep[]
    (mesh=out, matrix=[r.matrix[4 * (i - 1) + j] for i in 1:4, j in 1:4],
     mode=r.mode == 1 ? "axis" : "face", primary=collect(r.primary),
     secondary=collect(r.secondary), share=r.share, face_share=r.face_share,
     axis_share=r.axis_share, opening_area=r.opening_area)
end

const I4 = [1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0]

"""
    estimate_normals(m; k=16, orientation="consistent", viewpoint=(0,0,0),
                     record_variation=false) -> (; mesh, num_degenerate, num_components)

Normals of the mesh's points by local plane fits, as `normals` point data:
oriented `"consistent"` (a spanning tree over the neighbour graph, then each
piece outward), `"viewpoint"` (towards `viewpoint`) or `"none"`; `NaN` where a
neighbourhood spans no plane. See `doc/point_normals.md`.
"""
function estimate_normals(m::Mesh; k::Integer=16, orientation::AbstractString="consistent",
                          viewpoint=(0.0, 0.0, 0.0), record_variation::Bool=false)
    o = orientation == "consistent" ? Int32(0) : orientation == "viewpoint" ? Int32(1) :
        orientation == "none" ? Int32(2) :
        throw(ArgumentError("meshio++: estimate_normals: orientation must be \"consistent\", \"viewpoint\" or \"none\""))
    opts = _CPointNormalsOpts(Int64(k), o, record_variation ? Int32(1) : Int32(0),
                              (Float64(viewpoint[1]), Float64(viewpoint[2]), Float64(viewpoint[3])),
                              ntuple(_ -> Int64(0), 4))
    nd, nc = Ref{Int64}(0), Ref{Int64}(0)
    ptr = ccall(_sym(:mio_estimate_normals), Ptr{Cvoid},
                (Ptr{Cvoid}, Ref{_CPointNormalsOpts}, Ptr{Int64}, Ptr{Int64}),
                _handle(m), Ref(opts), nd, nc)
    r = _check_ptr(ptr)
    (mesh=Mesh(r), num_degenerate=Int(nd[]), num_components=Int(nc[]))
end

"""
    reconstruct_surface(m; normals="normals", k=8, cell_size=0.0, resolution=0,
                        max_distance=0.0, estimate_normals=true, normals_k=16,
                        max_cells=20_000_000)
        -> (; mesh, num_samples, spacing, cell_size, max_distance, num_lattice_points,
              num_undefined, num_triangles, num_components, estimated_normals)

A triangle surface through the mesh's points, an oriented point cloud: the zero
level set of the signed distance to the samples' tangent planes (Hoppe et al.
1992), cut by marching tetrahedra on a band of lattice cells about the samples.
The normals are the `normals` point data, estimated when absent unless
`estimate_normals=false`. See `doc/reconstruct_surface.md`.
"""
function reconstruct_surface(m::Mesh; normals::AbstractString="normals", k::Integer=8,
                             cell_size::Real=0.0, resolution::Integer=0, max_distance::Real=0.0,
                             estimate_normals::Bool=true, normals_k::Integer=16,
                             max_cells::Integer=20_000_000)
    rep = Ref{_CReconstructSurfaceReport}()
    ptr = GC.@preserve normals begin
        opts = _CReconstructSurfaceOpts(Base.unsafe_convert(Cstring, normals), Int64(k),
                                        Float64(cell_size), Int64(resolution),
                                        Float64(max_distance),
                                        estimate_normals ? Int32(1) : Int32(0), Int32(0),
                                        Int64(normals_k), Int64(max_cells),
                                        ntuple(_ -> Int64(0), 4))
        ccall(_sym(:mio_reconstruct_surface), Ptr{Cvoid},
              (Ptr{Cvoid}, Ref{_CReconstructSurfaceOpts}, Ref{_CReconstructSurfaceReport}),
              _handle(m), Ref(opts), rep)
    end
    r = _check_ptr(ptr)
    h = rep[]
    (mesh=Mesh(r), num_samples=Int(h.num_samples), spacing=h.spacing, cell_size=h.cell_size,
     max_distance=h.max_distance, num_lattice_points=Int(h.num_lattice_points),
     num_undefined=Int(h.num_undefined), num_triangles=Int(h.num_triangles),
     num_components=Int(h.num_components), estimated_normals=h.estimated_normals != 0)
end

"""
    convex_hull(m; dimension=0)
        -> (; mesh, dimension, num_vertices, num_facets, area, volume, perimeter)

The convex hull of the mesh's points (its cells are ignored): a closed outward
`triangle` surface in 3-D, a counter-clockwise `line` loop of the x, y
coordinates in 2-D. Only the extreme points are vertices, each with
`hull:point_id` (its 0-based input point). `dimension` is 2, 3, or 0 to decide
from the points. See `doc/hull.md`.
"""
function convex_hull(m::Mesh; dimension::Integer=0)
    rep = Ref{_CHullReport}()
    ptr = ccall(_sym(:mio_convex_hull), Ptr{Cvoid}, (Ptr{Cvoid}, Int32, Ref{_CHullReport}),
                _handle(m), Int32(dimension), rep)
    r = _check_ptr(ptr)
    h = rep[]
    (mesh=Mesh(r), dimension=Int(h.dimension), num_vertices=Int(h.num_vertices),
     num_facets=Int(h.num_facets), area=h.area, volume=h.volume, perimeter=h.perimeter)
end

"""
    delaunay_2d(m) -> (; mesh, num_triangles, num_duplicates, num_boundary)

The Delaunay triangulation of the mesh's points in the x, y plane: every point
kept (with its z and point data), one counter-clockwise `triangle` block; a
point coinciding in x, y with an earlier one is in no triangle. See
`doc/hull.md`.
"""
function delaunay_2d(m::Mesh)
    nt, nd, nb = Ref{Int64}(0), Ref{Int64}(0), Ref{Int64}(0)
    ptr = ccall(_sym(:mio_delaunay_2d), Ptr{Cvoid},
                (Ptr{Cvoid}, Ptr{Int64}, Ptr{Int64}, Ptr{Int64}), _handle(m), nt, nd, nb)
    r = _check_ptr(ptr)
    (mesh=Mesh(r), num_triangles=Int(nt[]), num_duplicates=Int(nd[]), num_boundary=Int(nb[]))
end

"""
    convert_coordinates(m; system="cylindrical", axis=(0,0,1), origin=(0,0,0),
                        reference=(1,0,0), arrays=String[], inverse=false,
                        record_coordinates=true)
        -> (; mesh, num_arrays, num_singular)

Cylindrical or spherical components about an axis: the points stay put, their
coordinates are recorded as point data (`coords:cylindrical` = `(r, theta, z)`,
`coords:spherical` = `(r, theta, phi)`), and each named vector or tensor array
(point or cell data) is rewritten in the local basis of its point, or back to
Cartesian with `inverse=true`. See `doc/coordinates.md`.
"""
function convert_coordinates(m::Mesh; system::AbstractString="cylindrical", axis=(0.0, 0.0, 1.0),
                             origin=(0.0, 0.0, 0.0), reference=(1.0, 0.0, 0.0),
                             arrays::AbstractVector{<:AbstractString}=String[],
                             inverse::Bool=false, record_coordinates::Bool=true)
    sys = system == "cylindrical" ? Int32(0) : system == "spherical" ? Int32(1) :
          throw(ArgumentError("meshio++: convert_coordinates: system must be \"cylindrical\" or \"spherical\""))
    names = [Vector{UInt8}(codeunits(String(a) * "\0")) for a in arrays]
    ptrs = [Cstring(pointer(n)) for n in names]
    na, ns = Ref{Int64}(0), Ref{Int64}(0)
    t3(v) = (Float64(v[1]), Float64(v[2]), Float64(v[3]))
    ptr = GC.@preserve names ptrs begin
        opts = _CCoordinateOpts(sys, inverse ? Int32(1) : Int32(0), t3(axis), t3(origin),
                                t3(reference),
                                isempty(ptrs) ? Ptr{Cstring}(C_NULL) : pointer(ptrs),
                                Int64(length(ptrs)), record_coordinates ? Int32(1) : Int32(0),
                                Int32(0), ntuple(_ -> Int64(0), 4))
        ccall(_sym(:mio_convert_coordinates), Ptr{Cvoid},
              (Ptr{Cvoid}, Ref{_CCoordinateOpts}, Ptr{Int64}, Ptr{Int64}),
              _handle(m), Ref(opts), na, ns)
    end
    r = _check_ptr(ptr)
    (mesh=Mesh(r), num_arrays=Int(na[]), num_singular=Int(ns[]))
end

"""
    geodesic_distance(m, sources=Int[]; source_region="", time_factor=1.0,
                      tolerance=1e-10, max_iterations=0)
        -> (; mesh, num_sources, num_unreached, num_degenerate, max_distance,
             time_step, heat_iterations, poisson_iterations, converged)

Geodesic distance on a triangle surface by the heat method: `mesh` carries the
point data `geodesic:distance`, the distance along the surface to the nearest
source. `sources` are **1-based** point ids; `source_region` names a region
whose points are sources too. A point in no triangle, or in a piece with no
source, gets `NaN`. See `doc/geodesic.md`.
"""
function geodesic_distance(m::Mesh, sources::AbstractVector{<:Integer}=Int[];
                           source_region::AbstractString="", time_factor::Real=1.0,
                           tolerance::Real=1e-10, max_iterations::Integer=0)
    report = Ref{_CGeodesicReport}()
    ids = Int64[Int64(s) - 1 for s in sources]
    region_c = Vector{UInt8}(codeunits(String(source_region) * "\0"))
    ptr = GC.@preserve ids region_c begin
        opts = _CGeodesicOpts(isempty(ids) ? Ptr{Int64}(C_NULL) : pointer(ids),
                              Int64(length(ids)), Cstring(pointer(region_c)),
                              Float64(time_factor), Float64(tolerance), Int64(max_iterations),
                              ntuple(_ -> Int64(0), 6))
        ccall(_sym(:mio_geodesic_distance), Ptr{Cvoid},
              (Ptr{Cvoid}, Ref{_CGeodesicOpts}, Ptr{_CGeodesicReport}),
              _handle(m), Ref(opts), report)
    end
    r = _check_ptr(ptr)
    rep = report[]
    (mesh=Mesh(r), num_sources=Int(rep.num_sources), num_unreached=Int(rep.num_unreached),
     num_degenerate=Int(rep.num_degenerate), max_distance=Float64(rep.max_distance),
     time_step=Float64(rep.time_step), heat_iterations=Int(rep.heat_iterations),
     poisson_iterations=Int(rep.poisson_iterations), converged=rep.converged != 0)
end
geodesic_distance(m::Mesh, source::Integer; kwargs...) = geodesic_distance(m, [source]; kwargs...)

"""
    parameterize(m; padding=0.02, max_depth=8, max_stretch=2.0)
        -> (; mesh, num_charts, num_flipped, stretch, packing_efficiency, num_seam_points)

A UV atlas of a triangle surface: `mesh` carries the point data `uv` (two
texture coordinates per point in the unit square), split at the seams
(`parameterize:parent_point` names each point's input point, 0-based), and the
cell data `parameterize:chart`. A chart that folds, or whose L2 stretch exceeds
`max_stretch`, is halved up to `max_depth` times. See `doc/parameterize.md`.
"""
function parameterize(m::Mesh; padding::Real=0.02, max_depth::Integer=8,
                      max_stretch::Real=2.0)
    report = Ref{_CParameterizeReport}()
    opts = _CParameterizeOpts(Float64(padding), Int64(max_depth), Float64(max_stretch),
                              ntuple(_ -> Int64(0), 5))
    ptr = ccall(_sym(:mio_parameterize), Ptr{Cvoid},
                (Ptr{Cvoid}, Ref{_CParameterizeOpts}, Ptr{_CParameterizeReport}),
                _handle(m), Ref(opts), report)
    r = _check_ptr(ptr)
    rep = report[]
    (mesh=Mesh(r), num_charts=Int(rep.num_charts), num_flipped=Int(rep.num_flipped),
     stretch=Float64(rep.stretch), packing_efficiency=Float64(rep.packing_efficiency),
     num_seam_points=Int(rep.num_seam_points))
end

"""
    surface_watertight_check(m) -> NamedTuple

What is wrong with a surface, in numbers: boundary edges, non-manifold edges,
inconsistently wound pairs and degenerate triangles. A signed distance is only
meaningful where these are zero.
"""
function surface_watertight_check(m::Mesh)
    out = Ref{_CSurfaceQuality}()
    _check(ccall(_sym(:mio_surface_watertight_check), Cint,
                 (Ptr{Cvoid}, Ptr{_CSurfaceQuality}), _handle(m), out))
    q = out[]
    (boundary_edges=Int(q.boundary_edges), non_manifold_edges=Int(q.non_manifold_edges),
     inconsistent_pairs=Int(q.inconsistent_pairs),
     degenerate_triangles=Int(q.degenerate_triangles), watertight=q.watertight != 0)
end

function _sdf_opts(; sign, weight, location, band, record_closest_cell, record_inside,
                   watertight_check, surface_region, grid_cell_size, max_winding_work,
                   region_buf)
    _CSdfOpts(Cstring(pointer(region_buf)), Cdouble(band), Cdouble(grid_cell_size),
              Cdouble(max_winding_work), _code(_SDF_SIGNS, sign, "sign"),
              _code(_SDF_WEIGHTS, weight, "weight"),
              _code(_SDF_LOCATIONS, location, "location"),
              _code(_SDF_CHECKS, watertight_check, "watertight check"),
              record_closest_cell ? Int32(1) : Int32(0),
              record_inside ? Int32(1) : Int32(0),
              (Int64(0), Int64(0), Int64(0), Int64(0), Int64(0), Int64(0)))
end

"""
    sample_distance(surface, points; sign=:pseudonormal, band=0.0, ...) -> Vector{Float64}

Signed distances from `points` -- a `3 x n` matrix, the column-major memory the
C ABI reads as `n` rows of three -- to `surface`. Negative is inside.
"""
function sample_distance(surface::Mesh, points::AbstractMatrix{<:Real};
                         sign=:pseudonormal, weight=:angle, band::Real=0.0,
                         watertight_check=:warn, surface_region::AbstractString="",
                         grid_cell_size::Real=0.0, max_winding_work::Real=2.0e9)
    size(points, 1) == 3 || throw(ArgumentError("points must be a 3 x n matrix"))
    n = size(points, 2)
    pts = Matrix{Cdouble}(points)
    out = Vector{Cdouble}(undef, max(n, 1))
    region_buf = Vector{UInt8}(codeunits(String(surface_region)))
    push!(region_buf, 0x00)
    GC.@preserve pts out region_buf begin
        opts = _sdf_opts(; sign, weight, location=:corner, band,
                         record_closest_cell=false, record_inside=false, watertight_check,
                         surface_region, grid_cell_size, max_winding_work, region_buf)
        _check(ccall(_sym(:mio_sample_distance), Cint,
                     (Ptr{Cvoid}, Ptr{Cdouble}, Int64, Ref{_CSdfOpts}, Ptr{Cdouble}),
                     _handle(surface), pointer(pts), Int64(n), Ref(opts), pointer(out)))
    end
    Float64[out[i] for i in 1:n]
end

"""
    distance_to_surface(query, surface; sign=:pseudonormal, location=:corner, ...) -> NamedTuple

Attach the signed distance from `query`'s points (or cell centres) to `surface`
as the `sdf:distance` array. Returns `(; mesh, num_banded, quality)`.
"""
function distance_to_surface(query::Mesh, surface::Mesh; sign=:pseudonormal, weight=:angle,
                             location=:corner, band::Real=0.0,
                             record_closest_cell::Bool=false, record_inside::Bool=false,
                             watertight_check=:warn, surface_region::AbstractString="",
                             grid_cell_size::Real=0.0, max_winding_work::Real=2.0e9)
    banded = Ref{Int64}(0)
    quality = Ref{_CSurfaceQuality}()
    region_buf = Vector{UInt8}(codeunits(String(surface_region)))
    push!(region_buf, 0x00)
    ptr = GC.@preserve region_buf begin
        opts = _sdf_opts(; sign, weight, location, band, record_closest_cell, record_inside,
                         watertight_check, surface_region, grid_cell_size, max_winding_work,
                         region_buf)
        ccall(_sym(:mio_distance_to_surface), Ptr{Cvoid},
              (Ptr{Cvoid}, Ptr{Cvoid}, Ref{_CSdfOpts}, Ptr{Int64}, Ptr{_CSurfaceQuality}),
              _handle(query), _handle(surface), Ref(opts), banded, quality)
    end
    q = quality[]
    (mesh=Mesh(_check_ptr(ptr)), num_banded=Int(banded[]),
     quality=(boundary_edges=Int(q.boundary_edges),
              non_manifold_edges=Int(q.non_manifold_edges),
              inconsistent_pairs=Int(q.inconsistent_pairs),
              degenerate_triangles=Int(q.degenerate_triangles),
              watertight=q.watertight != 0))
end

const _SDF_STRUCTURES = Dict("voxel" => Int32(0), "octree" => Int32(1))

"""
    compute_sdf(surface; structure=:voxel, resolution=nothing, cell_size=nothing,
                bounds=nothing, padding=0.0, padding_relative=0.1,
                root_resolution=8, max_depth=4, band_cells=1.0, record_levels=true,
                max_cells=20_000_000, sign=:pseudonormal, ...) -> NamedTuple

Generate a grid over `surface` and fill it with signed distances -- the one call
that turns a surface into a field.

`structure` is `:voxel` (a dense lattice) or `:octree` (refined near the
surface, and therefore **1-irregular**: it has hanging nodes; see
[`refine`](@ref)'s `closure=:balanced`). `resolution`/`cell_size` size a voxel
grid and are an **error** with `:octree`, whose finest cell is
`root_resolution / 2^max_depth` and is therefore already determined.

Returns `(; mesh, dims, origin, spacing, max_depth, num_banded, quality)`, where
`dims` are the ROOT cell counts and `spacing` the FINEST cell size. The mesh
also carries the numeric `sdf:*` `field_data` header; no file format persists
arbitrary `field_data`, so write it as `.vti`, whose Origin/Spacing/WholeExtent
attributes are the same information.
"""
function compute_sdf(surface::Mesh; structure=:voxel, resolution=nothing,
                     cell_size=nothing, bounds=nothing, padding::Real=0.0,
                     padding_relative::Real=0.1, root_resolution::Integer=8,
                     max_depth::Integer=4, band_cells::Real=1.0,
                     record_levels::Bool=true, max_cells::Integer=20_000_000,
                     sign=:pseudonormal, weight=:angle, location=:corner, band::Real=0.0,
                     record_closest_cell::Bool=false, record_inside::Bool=false,
                     watertight_check=:warn, surface_region::AbstractString="",
                     grid_cell_size::Real=0.0, max_winding_work::Real=2.0e9)
    res = resolution === nothing ? Int64[] : Int64[Int64(v) for v in resolution]
    (resolution === nothing || length(res) == 3) ||
        throw(ArgumentError("resolution must have three entries"))
    bnd = bounds === nothing ? Cdouble[] : Cdouble[Cdouble(v) for v in bounds]
    (bounds === nothing || length(bnd) == 6) ||
        throw(ArgumentError("bounds must have six entries"))

    dims = Vector{Int64}(undef, 3)
    origin = Vector{Cdouble}(undef, 3)
    spacing = Vector{Cdouble}(undef, 3)
    depth = Ref{Int64}(0)
    banded = Ref{Int64}(0)
    quality = Ref{_CSurfaceQuality}()
    region_buf = Vector{UInt8}(codeunits(String(surface_region)))
    push!(region_buf, 0x00)
    ptr = GC.@preserve res bnd dims origin spacing region_buf begin
        inner = _sdf_opts(; sign, weight, location, band, record_closest_cell,
                          record_inside, watertight_check, surface_region, grid_cell_size,
                          max_winding_work, region_buf)
        opts = _CComputeSdfOpts(isempty(res) ? Ptr{Int64}(C_NULL) : pointer(res),
                                isempty(bnd) ? Ptr{Cdouble}(C_NULL) : pointer(bnd),
                                Cdouble(cell_size === nothing ? 0.0 : cell_size),
                                Cdouble(padding), Cdouble(padding_relative),
                                Cdouble(band_cells), Int64(max_cells),
                                Int64(root_resolution), Int64(max_depth),
                                _code(_SDF_STRUCTURES, structure, "structure"),
                                record_levels ? Int32(1) : Int32(0),
                                (Int64(0), Int64(0), Int64(0), Int64(0), Int64(0), Int64(0)),
                                inner)
        ccall(_sym(:mio_compute_sdf), Ptr{Cvoid},
              (Ptr{Cvoid}, Ref{_CComputeSdfOpts}, Ptr{Int64}, Ptr{Cdouble}, Ptr{Cdouble},
               Ptr{Int64}, Ptr{Int64}, Ptr{_CSurfaceQuality}),
              _handle(surface), Ref(opts), pointer(dims), pointer(origin),
              pointer(spacing), depth, banded, quality)
    end
    q = quality[]
    (mesh=Mesh(_check_ptr(ptr)), dims=Tuple(Int.(dims)), origin=Tuple(Float64.(origin)),
     spacing=Tuple(Float64.(spacing)), max_depth=Int(depth[]), num_banded=Int(banded[]),
     quality=(boundary_edges=Int(q.boundary_edges),
              non_manifold_edges=Int(q.non_manifold_edges),
              inconsistent_pairs=Int(q.inconsistent_pairs),
              degenerate_triangles=Int(q.degenerate_triangles),
              watertight=q.watertight != 0))
end
