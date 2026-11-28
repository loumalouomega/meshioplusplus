"""
    SurfaceQuery(mesh; region="")

A reusable surface query handle (roadmap 9.3.2): the triangle soup and
hierarchies of `mesh`'s triangle, quad and polygon cells (or of cell `region`),
built once. [`cast_rays(::SurfaceQuery, ...)`](@ref) answers a batch of rays
exactly as `cast_rays(mesh, ...)` does, and [`closest_points`](@ref) the nearest
surface point of each query point. The handle copies the surface; it is
released by a finalizer, or deterministically by `close` (idempotent).
"""
mutable struct SurfaceQuery
    ptr::Ptr{Cvoid}
    function SurfaceQuery(ptr::Ptr{Cvoid})
        q = new(ptr)
        finalizer(q) do x
            if x.ptr != C_NULL
                ccall(_sym(:mio_surface_query_free), Cvoid, (Ptr{Cvoid},), x.ptr)
                x.ptr = C_NULL
            end
        end
        q
    end
end

SurfaceQuery(m::Mesh; region::AbstractString="") = SurfaceQuery(_check_ptr(ccall(
    _sym(:mio_surface_query_create), Ptr{Cvoid}, (Ptr{Cvoid}, Cstring), _handle(m), region)))

function _handle(q::SurfaceQuery)
    q.ptr == C_NULL && throw(MeshioError(MIO_ERR_INVALID_ARG, "surface query has been closed"))
    q.ptr
end

function Base.close(q::SurfaceQuery)
    if q.ptr != C_NULL
        ccall(_sym(:mio_surface_query_free), Cvoid, (Ptr{Cvoid},), q.ptr)
        q.ptr = C_NULL
    end
    nothing
end

"""
    num_triangles(q::SurfaceQuery) -> Int

The surface's triangles (quads and polygons fanned).
"""
num_triangles(q::SurfaceQuery) = Int(_check_count(ccall(_sym(:mio_surface_query_num_triangles),
                                                        Int64, (Ptr{Cvoid},), _handle(q)),
                                                  "mio_surface_query_num_triangles"))

"""
    cast_rays(q::SurfaceQuery, origins, directions; max_distance=0.0)
        -> (; hit, distance, point, normal, cell, num_hits)

`cast_rays(mesh, ...)` against the query's surface, bit-identically.
"""
function cast_rays(q::SurfaceQuery, origins::AbstractMatrix{<:Real},
                   directions::AbstractVecOrMat{<:Real}; max_distance::Real=0.0)
    size(origins, 1) == 3 || throw(ArgumentError("origins must be a 3 x n matrix"))
    dirs = directions isa AbstractVector ? reshape(Vector{Cdouble}(directions), 3, 1) :
           Matrix{Cdouble}(directions)
    size(dirs, 1) == 3 || throw(ArgumentError("directions must have three rows"))
    n = size(origins, 2)
    nd = size(dirs, 2)
    o = Matrix{Cdouble}(origins)
    hit = Vector{Int32}(undef, max(n, 1))
    distance = Vector{Cdouble}(undef, max(n, 1))
    point = Matrix{Cdouble}(undef, 3, max(n, 1))
    normal = Matrix{Cdouble}(undef, 3, max(n, 1))
    cell = Vector{Int64}(undef, max(n, 1))
    nh = Ref{Int64}(0)
    GC.@preserve o dirs hit distance point normal cell begin
        _check(ccall(_sym(:mio_surface_query_cast_rays), Cint,
                     (Ptr{Cvoid}, Ptr{Cdouble}, Int64, Ptr{Cdouble}, Int64, Cdouble,
                      Ptr{Int32}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Int64},
                      Ptr{Int64}),
                     _handle(q), pointer(o), Int64(n), pointer(dirs), Int64(nd),
                     Float64(max_distance), pointer(hit), pointer(distance), pointer(point),
                     pointer(normal), pointer(cell), nh))
    end
    (hit=Bool[hit[i] != 0 for i in 1:n], distance=distance[1:n], point=point[:, 1:n],
     normal=normal[:, 1:n], cell=Int[cell[i] + 1 for i in 1:n], num_hits=Int(nh[]))
end

"""
    closest_points(q::SurfaceQuery, points) -> (; point, distance, cell)

The nearest surface point of each column of the `3 x n` matrix `points`:
`point` (`3 x n`), `distance` (unsigned) and `cell` (1-based global input cell).
"""
function closest_points(q::SurfaceQuery, points::AbstractMatrix{<:Real})
    size(points, 1) == 3 || throw(ArgumentError("points must be a 3 x n matrix"))
    n = size(points, 2)
    p = Matrix{Cdouble}(points)
    point = Matrix{Cdouble}(undef, 3, max(n, 1))
    distance = Vector{Cdouble}(undef, max(n, 1))
    cell = Vector{Int64}(undef, max(n, 1))
    GC.@preserve p point distance cell begin
        _check(ccall(_sym(:mio_surface_query_closest_points), Cint,
                     (Ptr{Cvoid}, Ptr{Cdouble}, Int64, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Int64}),
                     _handle(q), pointer(p), Int64(n), pointer(point), pointer(distance),
                     pointer(cell)))
    end
    (point=point[:, 1:n], distance=distance[1:n], cell=Int[cell[i] + 1 for i in 1:n])
end
