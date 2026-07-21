!  License:         MIT License
!                   meshio++ default license: LICENSE
!
!  Main authors:    Vicente Mataix Ferrandiz
!
! Test program for the meshio++ Fortran module (registered as the `fortran_api`
! ctest). Plain checks, no framework: exit code 0 iff every check passed.
! argv(1) is a path PREFIX for the files the test writes (e.g.
! "<builddir>/fortran_test_out" -> "<builddir>/fortran_test_out_mesh.vtu").
!
! The mesh is deliberately non-square everywhere (5 points x 3 dims, 2 cells x
! 4 nodes, 3-component vector data) with asymmetric coordinates, so a
! transposed array mapping or a missed 1-based shift cannot cancel out.
program test_fortran_api
    use, intrinsic :: iso_fortran_env, only: real64, int64, error_unit
    use meshioplusplus
    implicit none

    type(mio_mesh) :: m, r, c, s, q, ro
    integer(int64) :: qcells, qinv, qdeg, bw
    integer(int64), allocatable :: node_perm(:)
    logical :: perm_ok, eq
    integer :: verdict
    real(real64) :: points(3, 5), vec(3, 5)
    real(real64), allocatable :: rpoints(:, :), rdata(:), rvec(:, :)
    real(real64), pointer :: pview(:, :)
    integer(int64) :: conn(4, 2)
    integer(int64), allocatable :: rconn(:, :)
    character(:), allocatable :: prefix, vtu_path, vtk_path, msg, name
    integer :: fails, ierr, nargs, arglen, i, j

    fails = 0
    nargs = command_argument_count()
    if (nargs < 1) then
        write (error_unit, '(a)') 'usage: test_fortran_api <output-path-prefix>'
        error stop 2
    end if
    call get_command_argument(1, length=arglen)
    allocate (character(arglen) :: prefix)
    call get_command_argument(1, value=prefix)
    vtu_path = prefix//'_mesh.vtu'
    vtk_path = prefix//'_mesh.vtk'

    ! ---- module-level introspection ------------------------------------
    call check(len(mio_version()) > 0, 'mio_version() is non-empty')
    call check(len(mio_mesh_backend()) > 0, 'mio_mesh_backend() is non-empty')
    call check(mio_format_readable('vtu'), 'vtu is readable')
    call check(mio_format_writable('vtu'), 'vtu is writable')
    call check(.not. mio_format_writable('openfoam'), 'openfoam is read-only')
    call check(.not. mio_format_readable('nonexistent'), 'unknown format is not readable')

    ! ---- build a small tet mesh from arrays ----------------------------
    points = reshape([0.0_real64, 0.0_real64, 0.0_real64, &
                      1.1_real64, 0.2_real64, 0.3_real64, &
                      0.4_real64, 1.2_real64, 0.5_real64, &
                      0.6_real64, 0.7_real64, 1.3_real64, &
                      1.4_real64, 1.5_real64, 1.6_real64], [3, 5])
    conn = reshape([1_int64, 2_int64, 3_int64, 4_int64, &
                    2_int64, 3_int64, 4_int64, 5_int64], [4, 2])
    do i = 1, 3
        do j = 1, 5
            vec(i, j) = 10.0_real64*j + i
        end do
    end do

    call m%set_points(points)
    call m%add_cell_block('tetra', conn)
    call m%add_point_data('temperature', [1.0_real64, 2.0_real64, 3.0_real64, 4.0_real64, &
                                          5.0_real64])
    call m%add_point_data('velocity', vec)
    call m%add_cell_data('quality', [0.5_real64, 0.75_real64])

    call check(m%num_points() == 5_int64, 'num_points before write')
    call check(m%point_dim() == 3_int64, 'point_dim before write')
    call check(m%num_cell_blocks() == 1_int64, 'num_cell_blocks before write')

    ! ---- mesh operations ------------------------------------------------
    s = m%extract_surface()
    call check(s%num_cell_blocks() >= 1_int64, 'extract_surface produces cells')
    call check(s%cell_block_type(1) == 'triangle', 'extract_surface gives triangles')
    call s%free()

    q = m%attach_quality()
    call check(q%num_cell_data() >= 1_int64, 'attach_quality adds cell_data')
    call q%free()

    call m%quality_counts(num_cells=qcells, num_inverted=qinv, num_degenerate=qdeg)
    call check(qcells == 2_int64, 'quality_counts reports two cells')

    ! ---- reorder + bandwidth -------------------------------------------
    bw = m%compute_bandwidth()
    call check(bw >= 0_int64, 'compute_bandwidth is non-negative')
    ro = m%reorder('rcm', node_perm=node_perm)
    call check(ro%num_points() == 5_int64, 'reorder preserves point count')
    call check(ro%num_cell_blocks() == 1_int64, 'reorder preserves cell blocks')
    call check(allocated(node_perm), 'reorder returns a node permutation')
    call check(size(node_perm) == 5, 'node permutation has one entry per node')
    ! 1-based bijection over [1, 5]
    perm_ok = .true.
    do i = 1, 5
        if (count(node_perm == int(i, int64)) /= 1) perm_ok = .false.
    end do
    call check(perm_ok, 'node permutation is a 1-based bijection')
    call ro%free()

    ! ---- write, read back, verify everything ---------------------------
    call m%write(vtu_path)
    call r%read(vtu_path)

    call check(r%num_points() == 5_int64, 'num_points after round-trip')
    call check(r%point_dim() == 3_int64, 'point_dim after round-trip')
    call check(r%num_cell_blocks() == 1_int64, 'num_cell_blocks after round-trip')
    call check(r%cell_block_type(1) == 'tetra', 'cell block type after round-trip')
    call check(r%cell_block_num_cells(1) == 2_int64, 'cell block num_cells')
    call check(r%cell_block_nodes_per_cell(1) == 4_int64, 'cell block nodes_per_cell')
    call check(.not. r%cell_block_is_ragged(1), 'cell block is not ragged')

    ! ---- diff / equals -------------------------------------------------
    eq = m%equals(r, atol=1.0e-8_real64)
    call check(eq, 'equals: round-trip mesh equals the original within tolerance')
    verdict = m%diff(r, atol=1.0e-8_real64)
    call check(verdict <= 1, 'diff: round-trip verdict is identical or within tolerance')

    call r%get_points(rpoints)
    call check(size(rpoints, 1) == 3 .and. size(rpoints, 2) == 5, 'points shape (dim, n)')
    call check(maxval(abs(rpoints - points)) < 1.0e-12_real64, 'point coordinates round-trip')

    call r%points_ptr(pview)
    call check(associated(pview), 'points_ptr is associated')
    call check(maxval(abs(pview - points)) < 1.0e-12_real64, 'points_ptr zero-copy view')

    call r%get_cell_block(1, rconn)
    call check(size(rconn, 1) == 4 .and. size(rconn, 2) == 2, 'connectivity shape (npc, ncells)')
    call check(all(rconn == conn), '1-based connectivity round-trip')

    call check(r%num_point_data() == 2_int64, 'num_point_data after round-trip')
    name = r%point_data_name(1)  ! sorted: "temperature" < "velocity"
    call check(name == 'temperature', 'first sorted point_data name')
    call r%get_point_data('temperature', rdata)
    call check(size(rdata) == 5, 'scalar point_data size')
    call check(maxval(abs(rdata - [1.0_real64, 2.0_real64, 3.0_real64, 4.0_real64, &
                                   5.0_real64])) < 1.0e-12_real64, 'scalar point_data values')
    call r%get_point_data('velocity', rvec)
    call check(size(rvec, 1) == 3 .and. size(rvec, 2) == 5, 'vector point_data shape (comp, n)')
    call check(maxval(abs(rvec - vec)) < 1.0e-12_real64, 'vector point_data values')

    call check(r%num_cell_data() == 1_int64, 'num_cell_data after round-trip')
    call check(r%cell_data_num_blocks('quality') == 1_int64, 'cell_data num_blocks')
    call r%get_cell_data('quality', 1, rdata)
    call check(maxval(abs(rdata - [0.5_real64, 0.75_real64])) < 1.0e-12_real64, &
               'cell_data values')

    ! ---- convert + read the result -------------------------------------
    call mio_convert(vtu_path, vtk_path)
    call c%read(vtk_path)
    call check(c%num_points() == 5_int64, 'num_points after convert to vtk')
    call check(c%num_cell_blocks() == 1_int64, 'num_cell_blocks after convert to vtk')

    ! ---- second-format round-trip (gmsh) exercises another writer/reader ----
    ! VTU is XML/binary; gmsh is a wholly different text writer/reader, so this
    ! shakes out format-specific bugs the single VTU path can't.
    block
        type(mio_mesh) :: g
        real(real64), allocatable :: gpoints(:, :)
        integer(int64), allocatable :: gconn(:, :)
        character(:), allocatable :: gmsh_path

        gmsh_path = prefix//'_mesh.msh'
        call m%write(gmsh_path)
        call g%read(gmsh_path)
        call check(g%num_points() == 5_int64, 'gmsh: num_points round-trip')
        call check(g%num_cell_blocks() == 1_int64, 'gmsh: num_cell_blocks round-trip')
        call check(g%cell_block_type(1) == 'tetra', 'gmsh: cell block type round-trip')
        call g%get_points(gpoints)
        call check(maxval(abs(gpoints - points)) < 1.0e-12_real64, &
                   'gmsh: point coordinates round-trip')
        call g%get_cell_block(1, gconn)
        call check(all(gconn == conn), 'gmsh: 1-based connectivity round-trip')
        call g%free()
    end block

    ! ---- error paths ----------------------------------------------------
    ierr = 0
    call r%read(prefix//'_does_not_exist.vtu', stat=ierr, errmsg=msg)
    call check(ierr /= 0, 'reading a nonexistent file sets stat')
    call check(len(msg) > 0, 'reading a nonexistent file sets errmsg')
    call check(len(mio_error_message()) > 0, 'mio_error_message() is populated')
    call check(r%num_points() == 5_int64, 'failed read leaves the previous mesh intact')

    call m%write(prefix//'_mesh.nonsense_extension', stat=ierr, errmsg=msg)
    call check(ierr /= 0, 'unknown extension sets stat')

    ! An explicit but unknown format name is rejected via stat/errmsg too.
    ierr = 0
    msg = ''
    call m%write(vtu_path, format='not-a-real-format', stat=ierr, errmsg=msg)
    call check(ierr /= 0, 'unknown explicit format sets stat')
    call check(len(msg) > 0, 'unknown explicit format sets errmsg')

    ! ---- data operations -----------------------------------------------
    ! These act on the data arrays only; the geometry must come through
    ! untouched. `m` carries point_data temperature/velocity and cell_data
    ! quality.
    block
        type(mio_mesh) :: d
        type(mio_data_array_info), allocatable :: arrays(:)
        character(len=STRBUF_LEN), allocatable :: dkeys(:)
        integer :: k
        logical :: saw

        ! drop / keep / rename
        d = m%data_drop(MIO_DATA_POINT, ['temperature'], stat=ierr)
        call check(ierr == 0, 'data_drop succeeds')
        call check(d%num_points() == 5_int64, 'data_drop leaves geometry intact')
        call check(d%num_point_data() == 1_int64, 'data_drop removed one array')
        call d%free()

        d = m%data_keep(MIO_DATA_POINT, ['velocity   '], stat=ierr)
        call check(ierr == 0, 'data_keep succeeds')
        call check(d%num_point_data() == 1_int64, 'data_keep retained one array')
        call d%free()

        d = m%data_rename(MIO_DATA_POINT, 'temperature', 'T', stat=ierr)
        call check(ierr == 0, 'data_rename succeeds')
        call check(d%num_point_data() == 2_int64, 'data_rename preserves the count')
        call d%free()

        ! An unknown key fails through stat rather than aborting.
        ierr = 0
        d = m%data_drop(MIO_DATA_POINT, ['nope'], stat=ierr, errmsg=msg)
        call check(ierr /= 0, 'data_drop on an unknown key sets stat')
        call check(len(msg) > 0, 'data_drop on an unknown key sets errmsg')

        ! averaging
        d = m%data_point_to_cell(stat=ierr)
        call check(ierr == 0, 'data_point_to_cell succeeds')
        call check(d%num_points() == 5_int64, 'data_point_to_cell keeps geometry')
        call d%free()

        d = m%data_cell_to_point(['quality'], weight=MIO_WEIGHT_MEASURE, stat=ierr)
        call check(ierr == 0, 'data_cell_to_point (measure-weighted) succeeds')
        call check(d%num_point_data() == 3_int64, 'data_cell_to_point added an array')
        call d%free()

        ! calc
        d = m%data_calc('norm(velocity)', 'speed', stat=ierr)
        call check(ierr == 0, 'data_calc succeeds')
        call check(d%num_point_data() == 3_int64, 'data_calc added an array')
        call d%free()

        ierr = 0
        d = m%data_calc('log(temperature)', 'bad', stat=ierr, errmsg=msg)
        call check(ierr /= 0, 'data_calc rejects an unknown function')
        call check(len(msg) > 0, 'data_calc failure sets errmsg')

        ! conditioning
        d = m%data_condition(MIO_DATA_POINT, ['temperature'], mode=MIO_COND_NORMALIZE, &
                             stat=ierr)
        call check(ierr == 0, 'data_condition normalize succeeds')
        call check(d%num_points() == 5_int64, 'data_condition keeps geometry')
        call d%free()

        ! read-only summary
        arrays = m%data_info(keys=dkeys, stat=ierr)
        call check(ierr == 0, 'data_info succeeds')
        call check(size(arrays) == 3, 'data_info reports every array')
        saw = .false.
        do k = 1, size(arrays)
            if (trim(dkeys(k)) == 'temperature') then
                saw = .true.
                call check(arrays(k)%location == MIO_DATA_POINT, 'temperature is point_data')
                call check(arrays(k)%num_entries == 5_int64, 'temperature has 5 entries')
                call check(arrays(k)%num_components == 1_int64, 'temperature is scalar')
                call check(arrays(k)%num_finite == 5_int64, 'temperature is all finite')
                call check(abs(arrays(k)%min - 1.0_real64) < 1.0e-12_real64, &
                           'temperature min is 1')
                call check(abs(arrays(k)%max - 5.0_real64) < 1.0e-12_real64, &
                           'temperature max is 5')
                call check(abs(arrays(k)%mean - 3.0_real64) < 1.0e-12_real64, &
                           'temperature mean is 3')
            end if
        end do
        call check(saw, 'data_info named the temperature array')
    end block

    ! -- selective reads and read_metadata ---------------------------------
    block
        type(mio_mesh) :: sel
        type(mio_metadata) :: meta
        integer :: st

        call m%write('fortran_selective.vtu', stat=st)
        call check(st == 0, 'write for selective read')

        ! points_only keeps geometry, drops every data array.
        call sel%read('fortran_selective.vtu', points_only=.true., stat=st)
        call check(st == 0, 'selective read succeeded')
        call check(sel%num_points() > 0, 'points_only kept the points')
        call sel%free()

        ! A zero-sized `arrays` means "no arrays", not "every array".
        block
            character(len=16) :: wanted(0)
            call sel%read('fortran_selective.vtu', arrays=wanted, stat=st)
            call check(st == 0, 'empty-arrays read succeeded')
            call sel%free()
        end block

        meta = mio_read_metadata('fortran_selective.vtu', stat=st)
        call check(st == 0, 'read_metadata succeeded')
        call check(meta%num_points > 0, 'metadata reports points')
        call check(allocated(meta%cell_blocks), 'metadata reports cell blocks')
        call check(size(meta%cell_blocks) > 0, 'metadata has >= 1 cell block')
        call check(meta%num_cells > 0, 'metadata reports cells')
        ! vtu has a native metadata path, so this was genuinely cheap.
        call check(.not. meta%fell_back_to_full_read, 'vtu metadata is not a fallback')
    end block

    ! -- convert_cells: elevate then linearize returns the original topology --
    block
        type(mio_mesh) :: up, down
        integer(int64), allocatable :: pmap(:)
        integer :: st

        up = m%convert_cells('elevate', stat=st)
        call check(st == 0, 'convert_cells elevate succeeded')
        call check(up%num_points() > m%num_points(), 'elevate added mid-edge nodes')

        down = up%convert_cells('linearize', point_map=pmap, stat=st)
        call check(st == 0, 'convert_cells linearize succeeded')
        call check(down%num_points() == m%num_points(), 'linearize pruned the mid nodes')
        call check(allocated(pmap), 'convert_cells returned a point map')
        call check(size(pmap) == up%num_points(), 'point map covers the input points')

        call up%free()
        call down%free()

        ! An unknown mode fails through stat, never by aborting.
        up = m%convert_cells('nope', stat=st)
        call check(st /= 0, 'convert_cells rejects an unknown mode')
    end block

    ! -- smooth: relax coordinates, leaving topology and data intact --
    block
        type(mio_mesh) :: relaxed, bad
        integer(int64) :: moved, skipped, n0, nc0
        real(real64) :: max_disp
        integer :: st

        n0 = m%num_points()
        nc0 = m%cell_block_num_cells(1)

        relaxed = m%smooth('taubin', 5, nodes_moved=moved, max_displacement=max_disp, &
                           skipped_inversion=skipped, stat=st)
        call check(st == 0, 'smooth succeeded')
        call check(relaxed%num_points() == n0, 'smooth preserves the point count')
        call check(relaxed%cell_block_num_cells(1) == nc0, 'smooth preserves the cell count')
        call check(moved >= 0_int64, 'smooth reported a node-moved count')
        call check(max_disp >= 0.0_real64, 'smooth reported a max displacement')
        call check(skipped >= 0_int64, 'smooth reported a skipped-inversion count')

        ! The optional arguments really are optional.
        call relaxed%free()
        relaxed = m%smooth('laplacian', 3, stat=st)
        call check(st == 0, 'smooth works without the optional out-args')

        ! An unknown method name fails through stat rather than aborting.
        bad = m%smooth('bogus', 1, stat=st)
        call check(st /= 0, 'smooth rejects an unknown method')

        call relaxed%free()
    end block

    ! -- refine: uniform subdivision into same-type children --
    block
        type(mio_mesh) :: one, two, quadratic, bad
        integer(int64), allocatable :: pmap(:)
        integer(int64) :: n0
        integer :: st

        n0 = m%cell_block_num_cells(1)

        one = m%refine(1, point_map=pmap, stat=st)
        call check(st == 0, 'refine succeeded')
        call check(one%cell_block_num_cells(1) == n0*8, 'refine split each tetra into 8')
        call check(one%num_points() > m%num_points(), 'refine added mid-edge nodes')
        call check(allocated(pmap), 'refine returned a point map')
        call check(size(pmap) == m%num_points(), 'point map covers the input points')

        ! levels=2 matches applying one level twice.
        two = m%refine(2, stat=st)
        call check(st == 0, 'refine levels=2 succeeded')
        call check(two%cell_block_num_cells(1) == n0*64, &
                   'two levels split each tetra into 64')

        ! A cell type with no same-type subdivision fails through stat.
        quadratic = m%convert_cells('elevate', stat=st)
        call check(st == 0, 'built a tetra10 mesh to reject')
        bad = quadratic%refine(1, stat=st)
        call check(st /= 0, 'refine rejects higher-order cells')

        call one%free()
        call two%free()
        call quadratic%free()
    end block

    ! -- partition: N balanced pieces + flat labels --
    block
        type(mio_mesh), allocatable :: pieces(:)
        integer(int64), allocatable :: labels(:)
        integer(int64) :: total_in, total_out, i
        integer :: st, p

        total_in = 0
        do p = 1, int(m%num_cell_blocks())
            total_in = total_in + m%cell_block_num_cells(p)
        end do

        pieces = m%partition(2, method='sfc', stat=st)
        call check(st == 0, 'partition succeeded')
        call check(size(pieces) == 2, 'partition returned exactly nparts pieces')
        total_out = 0
        do p = 1, size(pieces)
            call check(pieces(p)%num_cell_blocks() == m%num_cell_blocks(), &
                       'piece keeps the input block structure')
            total_out = total_out + pieces(p)%cell_block_num_cells(1)
        end do
        call check(total_out == total_in, 'pieces cover every cell exactly once')

        labels = m%partition_labels(2, method='sfc', stat=st)
        call check(st == 0, 'partition_labels succeeded')
        call check(size(labels, kind=int64) == total_in, 'labels cover every cell')
        call check(all(labels >= 0) .and. all(labels < 2), 'labels lie in [0, nparts)')

        do p = 1, size(pieces)
            call pieces(p)%free()
        end do
        deallocate (pieces)

        ! ghost_layers is a v1 stub and must fail through stat, not abort.
        pieces = m%partition(2, ghost_layers=1, stat=st)
        call check(st /= 0, 'partition rejects ghost_layers /= 0')
        do i = 1, size(pieces)
            call pieces(i)%free()
        end do
    end block

    call m%free()
    call r%free()
    call c%free()
    call check(.not. m%is_valid(), 'handle invalid after free')

    if (fails /= 0) then
        write (error_unit, '(a,i0,a)') 'test_fortran_api: ', fails, ' check(s) FAILED'
        error stop 1
    end if
    write (*, '(a)') 'test_fortran_api: all checks passed'

contains

    subroutine check(ok, what)
        logical, intent(in) :: ok
        character(*), intent(in) :: what
        if (.not. ok) then
            fails = fails + 1
            write (error_unit, '(a)') 'FAIL: '//what
        end if
    end subroutine

end program test_fortran_api
