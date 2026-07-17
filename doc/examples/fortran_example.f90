! Minimal meshio++ Fortran consumer: build a tet mesh, write it, read it back.
program example
    use, intrinsic :: iso_fortran_env, only: real64, int64
    use meshioplusplus
    implicit none
    type(mio_mesh) :: m, r
    real(real64) :: points(3, 5)
    integer(int64) :: conn(4, 2)
    integer(int64), allocatable :: rconn(:, :)

    print '(a)', 'meshio++ '//mio_version()//' (backend: '//mio_mesh_backend()//')'

    points = reshape([0.0_real64, 0.0_real64, 0.0_real64, &
                      1.1_real64, 0.2_real64, 0.3_real64, &
                      0.4_real64, 1.2_real64, 0.5_real64, &
                      0.6_real64, 0.7_real64, 1.3_real64, &
                      1.4_real64, 1.5_real64, 1.6_real64], [3, 5])
    conn = reshape([1_int64, 2_int64, 3_int64, 4_int64, &
                    2_int64, 3_int64, 4_int64, 5_int64], [4, 2])

    call m%set_points(points)
    call m%add_cell_block('tetra', conn)
    call m%write('/tmp/mio_example_f.vtu')
    call m%free()

    call r%read('/tmp/mio_example_f.vtu')
    if (r%num_points() /= 5_int64) error stop 'wrong point count'
    call r%get_cell_block(1, rconn)
    if (.not. all(rconn == conn)) error stop 'wrong connectivity'
    call r%free()
    print '(a)', 'example.f90: OK'
end program example
