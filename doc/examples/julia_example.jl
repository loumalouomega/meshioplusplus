# Minimal meshio++ Julia consumer: build a tet mesh, write it, read it back.
#
# Run with the C library discoverable, e.g.
#   MESHIOPLUSPLUS_LIB=/opt/meshioplusplus/lib/libmeshioplusplus.so \
#     julia --project=bindings/julia/MeshioPlusPlus doc/examples/julia_example.jl

using MeshioPlusPlus
import MeshioPlusPlus as mio

println("meshio++ $(version()) (backend: $(mesh_backend()))")

# (dim, num_points) and (nodes_per_cell, num_cells): column-major Julia is the
# same memory as the C API's row-major (num_points, dim) / (num_cells, npc), so
# nothing is transposed. Connectivity is 1-based here.
pts = Float64[0.0 1.1 0.4 0.6 1.4
                 0.0 0.2 1.2 0.7 1.5
                 0.0 0.3 0.5 1.3 1.6]
conn = Int64[1 2
             2 3
             3 4
             4 5]

path = joinpath(tempdir(), "mio_example_julia.vtu")

m = Mesh()
set_points!(m, pts)
add_cell_block!(m, "tetra", conn)
mio.write(m, path)
close(m)

r = mio.read(path)
try
    num_points(r) == 5 || error("verification failed: point count")
    points(r) ≈ pts || error("verification failed: coordinates")
    cell_block_type(r, 1) == "tetra" || error("verification failed: cell type")
    connectivity(r, 1) == conn || error("verification failed: 1-based connectivity")
    # The borrow is 0-based -- the ABI's own numbering, deliberately unshifted.
    parent(connectivity_ptr(r, 1)) == conn .- 1 || error("verification failed: borrow")
finally
    close(r)
    rm(path; force=true)
end

println("julia_example.jl: OK")
