# Test suite for the meshio++ Julia binding.
#
# The fixture mirrors tests/fortran/test_fortran_api.f90 deliberately: 5 points
# x 3 dims, 2 tetra x 4 nodes, 3-component vector data, and coordinates that
# are asymmetric in every axis. Nothing about it is square, so a transposed
# array mapping or a missed 1-based shift cannot cancel out and pass anyway.

using Test
using MeshioPlusPlus
import MeshioPlusPlus as mio

# The same numbers the Fortran test uses, in Julia's (dim, num_points) shape.
const POINTS = Float64[0.0 1.1 0.4 0.6 1.4
                       0.0 0.2 1.2 0.7 1.5
                       0.0 0.3 0.5 1.3 1.6]
# (nodes_per_cell, num_cells), 1-based.
const CONN = Int64[1 2
                   2 3
                   3 4
                   4 5]
const VEC = Float64[10.0 20.0 30.0 40.0 50.0
                    11.0 21.0 31.0 41.0 51.0
                    12.0 22.0 32.0 42.0 52.0]

"""Build the reference mesh from raw arrays."""
function fixture()
    m = Mesh()
    set_points!(m, POINTS)
    add_cell_block!(m, "tetra", CONN)
    add_point_data!(m, "displacement", VEC)
    add_point_data!(m, "temperature", Float64[1.0, 2.0, 3.0, 4.0, 5.0])
    append_cell_data!(m, "material", Int64[7, 9])
    m
end

@testset "MeshioPlusPlus" begin

@testset "introspection" begin
    @test !isempty(version())
    @test mesh_backend() in ("meshio", "native", "kratos")
    @test !isempty(library_path())
    @test format_readable("vtu")
    @test format_writable("vtu")
    @test !format_writable("openfoam")      # read-only format
    @test !format_readable("nonexistent")
    @test cell_type_num_nodes("tetra10") == 10
    @test cell_type_dimension("triangle") == 2
    @test cell_type_num_nodes("polygon") == -1   # variable-size
end

@testset "build and inspect" begin
    m = fixture()
    @test num_points(m) == 5
    @test point_dim(m) == 3
    @test num_cell_blocks(m) == 1
    @test num_cells(m) == 2
    @test cell_block_type(m, 1) == "tetra"
    @test cell_block_types(m) == ["tetra"]

    info = cell_block_info(m, 1)
    @test info.num_cells == 2
    @test info.nodes_per_cell == 4
    @test !info.is_ragged

    # Data names come back sorted, on every backend.
    @test point_data_names(m) == ["displacement", "temperature"]
    @test cell_data_names(m) == ["material"]
    @test num_point_data(m) == 2
    @test cell_data_num_blocks(m, "material") == 1
    close(m)
end

@testset "column-major shape identity" begin
    m = fixture()
    p = points(m)
    # (dim, num_points) in Julia == the C API's row-major (num_points, dim):
    # the SAME memory, never transposed. If a transpose had crept in, this
    # non-square 3x5 fixture would come back 5x3 or scrambled.
    @test size(p) == (3, 5)
    @test p == POINTS
    @test p[:, 2] == [1.1, 0.2, 0.3]      # point 2's coordinates
    @test p[1, :] == [0.0, 1.1, 0.4, 0.6, 1.4]  # every point's x

    # A 3-component vector field is (components, num_points) here.
    d = point_data(m, "displacement")
    @test size(d) == (3, 5)
    @test d == VEC
    # A scalar field stays a plain vector.
    @test point_data(m, "temperature") == [1.0, 2.0, 3.0, 4.0, 5.0]
    @test cell_data(m, "material", 1) == Int64[7, 9]
    close(m)
end

@testset "1-based copying vs 0-based borrowing" begin
    m = fixture()
    c1 = connectivity(m, 1)             # copy, 1-based
    c0 = parent(connectivity_ptr(m, 1)) # borrow, 0-based (the ABI's own)
    @test size(c1) == (4, 2)
    @test c1 == CONN
    @test c1 == c0 .+ 1                 # exactly the documented relationship
    @test eltype(c0) == Int64
    close(m)
end

@testset "zero-copy borrow really aliases" begin
    m = fixture()
    b = points_ptr(m)
    @test size(b) == (3, 5)
    @test b[1, 2] == 1.1
    # Writing through the borrow must be visible to the mesh, or it is not a
    # borrow at all.
    b[1, 2] = 99.0
    @test points(m)[1, 2] == 99.0
    close(m)
end

@testset "borrow window is enforced" begin
    m = fixture()
    b = points_ptr(m)
    @test b[1, 2] == 1.1                     # readable inside the window

    # A MUTATING call ends the window. The C API says a borrowed pointer stays
    # valid only until the next mutating mio_mesh_* call; rather than let a
    # stale read succeed by luck, using the borrow afterwards must throw.
    add_point_data!(m, "extra", Float64[1.0, 2.0, 3.0, 4.0, 5.0])
    @test_throws BorrowError b[1, 2]
    @test_throws BorrowError size(b)

    # A fresh borrow works again.
    b2 = points_ptr(m)
    @test b2[1, 2] == 1.1

    # Closing the mesh also ends every window.
    close(m)
    @test_throws BorrowError b2[1, 2]

    # A read-only accessor does NOT end the window (C API rule 3).
    m2 = fixture()
    b3 = points_ptr(m2)
    _ = num_points(m2); _ = points(m2); _ = connectivity(m2, 1)
    @test b3[1, 2] == 1.1
    close(m2)
end

@testset "file round-trip" begin
    mktempdir() do dir
        m = fixture()
        for ext in ("vtu", "vtk")
            path = joinpath(dir, "mesh.$ext")
            mio.write(m, path)
            @test isfile(path)
            r = mio.read(path)
            @test num_points(r) == 5
            @test point_dim(r) == 3
            @test cell_block_type(r, 1) == "tetra"
            @test points(r) ≈ POINTS
            @test connectivity(r, 1) == CONN
            @test point_data(r, "displacement") ≈ VEC
            @test meshes_equal(m, r; atol=1e-10, rtol=1e-10)
            close(r)
        end

        # convert() never materializes a mesh for the caller.
        src = joinpath(dir, "mesh.vtu")
        dst = joinpath(dir, "converted.vtk")
        mio.convert(src, dst)
        @test isfile(dst)

        # sniff_format reads the leading bytes, not the extension.
        @test mio.sniff_format(src) == "vtu"
        close(m)
    end
end

@testset "selective read and metadata" begin
    mktempdir() do dir
        path = joinpath(dir, "mesh.vtu")
        m = fixture()
        mio.write(m, path)
        close(m)

        r = mio.read(path; options=ReadOptions(points_only=true))
        @test num_points(r) == 5
        @test isempty(point_data_names(r))    # every array skipped
        close(r)

        # `arrays=[...]` keeps only what is named; the empty vector is NOT the
        # same request as `nothing` (which keeps everything).
        r2 = mio.read(path; options=ReadOptions(arrays=["temperature"]))
        @test point_data_names(r2) == ["temperature"]
        close(r2)
        r3 = mio.read(path; options=ReadOptions(arrays=String[]))
        @test isempty(point_data_names(r3))
        close(r3)

        meta = read_metadata(path)
        @test meta.num_points == 5
        @test meta.point_dim == 3
        @test meta.num_cells == 2
        @test meta.num_cell_blocks == 1
        @test meta.cell_blocks[1].type == "tetra"
        @test meta.cell_blocks[1].num_cells == 2
        @test !meta.cell_blocks[1].is_ragged
        @test "temperature" in meta.point_data_names
        @test reader_supports_options("vtu")
    end
end

@testset "regions round-trip" begin
    m = fixture()
    add_region!(m, "inlet", :point, [1, 3, 5])
    add_region!(m, "solid", :cell, [1, 2]; dim=3, tag=17)
    # (cell, facet) pairs: column-major, so 2 rows x N columns.
    add_region!(m, "wall", :side, Int64[1 2; 0 2])

    rs = regions(m)
    @test length(rs) == 3
    by_name = Dict(r.name => r for r in rs)

    inlet = by_name["inlet"]
    @test inlet.kind === :point
    @test vec(inlet.entries) == [1, 3, 5]     # 1-based, as given

    solid = by_name["solid"]
    @test solid.kind === :cell
    @test solid.dim == 3
    @test solid.tag == 17
    @test vec(solid.entries) == [1, 2]

    wall = by_name["wall"]
    @test wall.kind === :side
    @test size(wall.entries) == (2, 2)
    # Canonicalized (sorted lexicographically on the pair) by the core. The
    # CELL column is 1-based here; the FACET column is an ordinal within the
    # cell type, not a mesh index, so it is NOT shifted -- the same rule the
    # Fortran module follows.
    @test wall.entries[1, :] == [1, 2]        # cells, 1-based
    @test wall.entries[2, :] == [0, 2]        # facets, unshifted

    @test_throws ArgumentError add_region!(m, "bad", :point, [0, 1])   # 0 is not 1-based
    @test_throws ArgumentError add_region!(m, "bad", :side, [1, 2])    # needs pairs
    close(m)
end

@testset "regions survive a file round-trip" begin
    mktempdir() do dir
        m = Mesh()
        set_points!(m, POINTS)
        add_cell_block!(m, "triangle", Int64[1 2; 2 3; 3 4])
        add_region!(m, "physical", :cell, [1, 2]; dim=2, tag=5)

        # Abaqus is one of the two Phase-1 region formats. Its *ELSET has no
        # notion of a dimension or an integer tag, so those come back
        # unspecified (-1) -- the group itself is what survives.
        path = joinpath(dir, "mesh.inp")
        mio.write(m, path)
        r = mio.read(path)
        got = only(filter(x -> x.name == "physical", regions(r)))
        @test got.kind === :cell
        @test vec(got.entries) == [1, 2]
        close(r); close(m)
    end
end

@testset "operations: geometry" begin
    m = fixture()

    surf = extract_surface(m)
    @test num_points(surf) > 0
    @test cell_block_type(surf, 1) == "triangle"
    close(surf)

    surf2 = extract_surface(m; record_parent_ids=true)
    @test "surface:parent_cell" in cell_data_names(surf2)
    close(surf2)

    skin = extract_skin(m; linearize=true)
    @test num_cells(skin) > 0
    close(skin)

    q = attach_quality(m)
    @test any(startswith(n, "quality:") for n in cell_data_names(q))
    close(q)

    qc = quality_counts(m)
    @test qc.num_cells == 2

    st = stats(m)
    @test st.num_points == 5
    @test st.num_cells == 2
    @test st.bbox_min[1] ≈ 0.0
    @test st.bbox_max[1] ≈ 1.4
    @test st.extent[1] ≈ 1.4

    @test compute_bandwidth(m) >= 0
    close(m)
end

@testset "operations: transform" begin
    m = fixture()
    M = [1.0 0 0 2.0
         0 1.0 0 3.0
         0 0 1.0 4.0
         0 0 0 1.0]
    t = mio.transform(m, M)
    p = points(t)
    @test p[:, 1] ≈ [2.0, 3.0, 4.0]          # the origin moved by the translation
    @test p[:, 2] ≈ POINTS[:, 2] .+ [2.0, 3.0, 4.0]
    @test connectivity(t, 1) == CONN         # topology untouched
    close(t); close(m)
    @test_throws ArgumentError mio.transform(fixture(), [1.0 0; 0 1.0])
end

@testset "operations: refine and the map sentinel" begin
    m = fixture()
    r = refine(m; levels=1)
    @test num_cells(r.mesh) == 16            # 2 tetra -> 8 children each
    # Refinement never prunes, so the point map is the identity -- and it is
    # 1-based here, like every copied index array.
    @test r.point_map == collect(1:5)
    @test length(r.cell_maps) == 1
    @test r.cell_maps[1] == [1, 9]           # each parent's first child, 1-based
    close(r.mesh)

    # convert_cells CAN prune, and a pruned point is reported as 0 -- never a
    # valid 1-based index, so the sentinel stays unambiguous (Fortran's rule).
    c = convert_cells(m, "simplexify")
    @test num_cells(c.mesh) == 2             # a tetra is already a simplex
    @test all(>=(0), c.point_map)
    close(c.mesh)
    close(m)
end

@testset "operations: selective refine with a conforming closure" begin
    m = fixture()

    # One cell only. `cells` is 1-based here, like every index array this
    # binding takes, and is shifted to the C API's 0-based numbering inside.
    sel = refine(m; cells=[1], record_levels=true)
    @test num_cells(sel.mesh) > num_cells(m)
    @test num_cells(sel.mesh) < 16           # not the uniform 8-per-cell
    @test "refine:level" in cell_data_names(sel.mesh)
    close(sel.mesh)

    # Propagation is the always-works baseline: on a connected mesh it reaches
    # every cell, which is exactly the uniform refinement.
    prop = refine(m; cells=[1], closure="propagate")
    @test num_cells(prop.mesh) == 16
    close(prop.mesh)

    # At most one selector, and the closure name is validated before the call.
    @test_throws Exception refine(m; cells=[1], region="anything")
    @test_throws ArgumentError refine(m; closure="blue")
    @test_throws ArgumentError refine(m; where_array="q", where_op="~=")
    close(m)
end

@testset "operations: split and partition own their meshes" begin
    m = fixture()

    pieces = mio.split(m; by="type")
    @test length(pieces) == 1
    key, piece = pieces[1]
    @test key == "tetra"
    # The piece must still be usable after the C result handle is gone --
    # which it is, because every piece takes ownership rather than borrowing.
    @test num_cells(piece) == 2
    @test points(piece) == POINTS
    close(piece)

    parts = partition(m, 2; method="sfc")
    @test length(parts) == 2
    @test [p.part_id for p in parts] == [0, 1]
    @test sum(num_cells(p.mesh) for p in parts) == 2   # partition of unity
    for p in parts
        @test length(p.cell_maps) == 1                 # block structure kept 1:1
        close(p.mesh)
    end

    labels = partition_labels(m, 2; method="sfc")
    @test length(labels) == 2
    # Part IDs, not indices: deliberately NOT shifted to 1-based.
    @test all(l -> 0 <= l < 2, labels)
    @test sort(labels) == [0, 1]
    close(m)
end

@testset "operations: reorder" begin
    m = fixture()
    r = reorder(m, "rcm")
    @test num_points(r.mesh) == 5
    @test sort(r.node_perm) == collect(1:5)   # a permutation, 1-based
    @test length(r.cell_perms) == 1
    @test sort(r.cell_perms[1]) == collect(1:2)
    close(r.mesh)
    close(m)
end

@testset "operations: clean, smooth, crop" begin
    m = fixture()

    c = clean(m)
    @test num_points(c.mesh) == 5             # nothing to remove in the fixture
    @test c.points_welded == 0
    close(c.mesh)

    s = smooth(m; iterations=2)
    @test num_points(s.mesh) == 5             # a pure coordinate move
    @test num_cells(s.mesh) == 2
    @test s.max_displacement >= 0.0
    close(s.mesh)

    b = crop_bbox(m, [-1.0, -1.0, -1.0], [10.0, 10.0, 10.0])
    @test num_cells(b) == 2                   # the box holds everything
    close(b)
    b2 = crop_bbox(m, [100.0, 100.0, 100.0], [101.0, 101.0, 101.0])
    @test num_cells(b2) == 0
    close(b2)

    p = crop_plane(m, [0.0, 0.0, 0.0], [1.0, 0.0, 0.0]; mode=:any)
    @test num_cells(p) >= 0
    close(p)
    close(m)
end

@testset "operations: slice and isosurface" begin
    m = fixture()

    sl = mio.slice(m, [0.5, 0.5, 0.5], [0.0, 0.0, 1.0])
    @test num_points(sl) >= 0                 # may be empty; never an error
    close(sl)

    iso = isosurface(m, "temperature", 3.0)
    @test num_points(iso) >= 0
    if num_cells(iso) > 0
        @test "iso:value" in cell_data_names(iso)
        @test "iso:index" in cell_data_names(iso)
    end
    close(iso)

    # A cell_data array has no level set: that must fail by name.
    @test_throws MeshioError isosurface(m, "material", 8.0)
    close(m)
end

@testset "operations: gradient" begin
    m = fixture()

    g = gradient(m, "temperature")
    @test g.num_skipped == 0
    @test "temperature:gradient" in cell_data_names(g.mesh)
    close(g.mesh)

    # The underscore spelling must reach the C API's hyphenated name.
    l = gradient(m, "temperature"; method=:least_squares)
    @test l.num_fallback >= 0
    close(l.mesh)

    p = gradient(m, "temperature"; location=:point, output="dT")
    @test "dT" in point_data_names(p.mesh)
    close(p.mesh)

    # A cell_data array is piecewise constant and has no derivative, and a
    # scalar has no divergence: both must fail by name.
    @test_throws MeshioError gradient(m, "material")
    @test_throws MeshioError gradient(m, "temperature"; operator=:divergence)
    close(m)
end

@testset "operations: merge and interpolate" begin
    a = fixture()
    b = fixture()
    merged = mio.merge([a, b])
    @test num_points(merged) == 10
    @test num_cells(merged) == 4
    @test "source_mesh_id" in cell_data_names(merged)
    close(merged)

    merged2 = mio.merge([a, b]; weld=true, atol=1e-9, source_tag=false)
    @test num_points(merged2) == 5            # the two copies weld together
    close(merged2)

    # The target already carries "temperature", so the default on_conflict
    # ("error") must refuse rather than silently overwrite.
    target = fixture()
    @test_throws MeshioError interpolate(a, target; arrays=["temperature"])

    it = interpolate(a, target; arrays=["temperature"], on_conflict="suffix")
    @test "temperature_interp" in point_data_names(it)
    @test point_data(it, "temperature_interp") ≈ [1.0, 2.0, 3.0, 4.0, 5.0]
    close(it)

    # Onto a target with no data of its own, the plain call works.
    bare = Mesh()
    set_points!(bare, POINTS)
    add_cell_block!(bare, "tetra", CONN)
    it2 = interpolate(a, bare; arrays=["temperature"])
    @test point_data(it2, "temperature") ≈ [1.0, 2.0, 3.0, 4.0, 5.0]
    close(it2); close(bare); close(target); close(a); close(b)
end

@testset "operations: diff" begin
    a = fixture()
    b = fixture()
    rep = mio.diff(a, b)
    @test rep.verdict === :identical
    @test rep.points.num_exceeding == 0

    pb = points_ptr(b)
    pb[1, 1] = 5.0
    rep2 = mio.diff(a, b)
    @test rep2.verdict === :different
    @test rep2.points.max_abs_error ≈ 5.0
    @test !meshes_equal(a, b)
    @test meshes_equal(a, b; atol=10.0)       # inside a loose tolerance
    close(a); close(b)
end

@testset "data operations never touch geometry" begin
    m = fixture()
    before_points = points(m)
    before_conn = connectivity(m, 1)

    d = data_drop(m, :point, ["temperature"])
    @test point_data_names(d) == ["displacement"]
    @test points(d) == before_points          # geometry bit-identical
    @test connectivity(d, 1) == before_conn
    close(d)

    k = data_keep(m, :point, ["temperature"])
    @test point_data_names(k) == ["temperature"]
    close(k)

    r = data_rename(m, :point, "temperature", "T")
    @test "T" in point_data_names(r)
    @test !("temperature" in point_data_names(r))
    close(r)

    p2c = data_point_to_cell(m, ["temperature"])
    @test "temperature" in cell_data_names(p2c)
    close(p2c)

    c2p = data_cell_to_point(m, ["material"]; weight=:uniform)
    @test "material" in point_data_names(c2p)
    close(c2p)

    calc = data_calc(m, "temperature * 2.0", "doubled")
    @test point_data(calc, "doubled") ≈ [2.0, 4.0, 6.0, 8.0, 10.0]
    close(calc)

    cond = data_condition(m, :point, ["temperature"]; mode=:normalize, lo=0.0, hi=1.0)
    t = point_data(cond, "temperature")
    @test minimum(t) ≈ 0.0
    @test maximum(t) ≈ 1.0
    close(cond)

    info = data_info(m)
    @test length(info) == 3
    temp = only(filter(x -> x.name == "temperature", info))
    @test temp.location === :point
    @test temp.num_entries == 5
    @test temp.num_components == 1
    @test temp.min ≈ 1.0
    @test temp.max ≈ 5.0
    @test temp.mean ≈ 3.0
    @test temp.num_nan == 0
    @test length(temp.components) == 1
    close(m)
end

@testset "transient XDMF series" begin
    # The one writer `write` cannot express: the grid goes out once and each
    # step is appended. "XML" keeps everything in the single .xdmf, so this
    # runs against a library built without HDF5 too.
    mktempdir() do dir
        path = joinpath(dir, "series.xdmf")
        m = fixture()

        s = XdmfSeries(path; data_format="XML")
        @test isopen(s)
        @test num_steps(s) == 0
        write_points_cells!(s, m)

        # temperature = t + node id, so no two steps share a value and a
        # step mix-up cannot pass.
        step = Mesh()
        set_points!(step, POINTS)
        add_cell_block!(step, "tetra", CONN)
        times = [0.0, 0.5, 1.0]
        for t in times
            add_point_data!(step, "temperature", Float64[t + i for i in 1:5])
            write_data!(s, t, step)
        end
        @test num_steps(s) == 3
        @test occursin("3 steps", sprint(show, s))

        # The .xdmf is buffered until finalize; nothing is readable before it.
        finalize!(s)
        finalize!(s)                       # idempotent
        close(s)
        @test !isopen(s)
        close(s)                           # idempotent
        @test_throws MeshioError num_steps(s)
        close(step)

        meta = read_metadata(path)
        @test meta.time_values ≈ times

        for (k, t) in enumerate(times)
            back = mio.read(path; options=ReadOptions(time_step=k - 1))
            @test num_points(back) == 5
            @test point_data(back, "temperature") ≈ Float64[t + i for i in 1:5]
            close(back)
        end

        # The do-block form closes even when the body throws.
        outer = nothing
        @test_throws ErrorException XdmfSeries(joinpath(dir, "boom.xdmf");
                                               data_format="XML") do inner
            outer = inner
            write_points_cells!(inner, m)
            error("boom")
        end
        @test !isopen(outer)

        # An unknown data format throws from the constructor, carrying the C
        # API's own message.
        @test_throws MeshioError XdmfSeries(joinpath(dir, "bad.xdmf");
                                            data_format="NoSuchFormat")
        close(m)
    end
end

@testset "errors carry mio_last_error()" begin
    # Every failure surfaces as a MeshioError carrying the C API's own
    # thread-local message; a status code never reaches the caller.
    err = try
        mio.read("/definitely/not/a/mesh/file.vtu")
        nothing
    catch e
        e
    end
    @test err isa MeshioError
    @test !isempty(err.msg)
    @test err.msg == last_error()             # straight from mio_last_error()
    @test occursin("MeshioError", sprint(showerror, err))

    m = fixture()
    @test_throws MeshioError point_data(m, "no_such_array")
    @test_throws MeshioError cell_block_type(m, 99)
    @test_throws MeshioError data_rename(m, :point, "missing", "x")

    # Using a closed mesh is an error, not a crash.
    close(m)
    @test_throws MeshioError num_points(m)
    close(m)                                   # idempotent
    @test !isopen(m)
end

end # testset
