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

    type(mio_mesh) :: m, r, c
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

    ! ---- error paths ----------------------------------------------------
    ierr = 0
    call r%read(prefix//'_does_not_exist.vtu', stat=ierr, errmsg=msg)
    call check(ierr /= 0, 'reading a nonexistent file sets stat')
    call check(len(msg) > 0, 'reading a nonexistent file sets errmsg')
    call check(len(mio_error_message()) > 0, 'mio_error_message() is populated')
    call check(r%num_points() == 5_int64, 'failed read leaves the previous mesh intact')

    call m%write(prefix//'_mesh.nonsense_extension', stat=ierr, errmsg=msg)
    call check(ierr /= 0, 'unknown extension sets stat')

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
