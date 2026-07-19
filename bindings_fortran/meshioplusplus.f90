!  ██████   ██████ ██████████  █████████  █████   █████ █████    ███████
! ░░██████ ██████ ░░███░░░░░█ ███░░░░░███░░███   ░░███ ░░███   ███░░░░░███      ███         ███
!  ░███░█████░███  ░███  █ ░ ░███    ░░░  ░███    ░███  ░███  ███     ░░███    ░███        ░███
!  ░███░░███ ░███  ░██████   ░░█████████  ░███████████  ░███ ░███      ░███ ███████████ ███████████
!  ░███ ░░░  ░███  ░███░░█    ░░░░░░░░███ ░███░░░░░███  ░███ ░███      ░███░░░░░███░░░ ░░░░░███░░░
!  ░███      ░███  ░███ ░   █ ███    ░███ ░███    ░███  ░███ ░░███     ███     ░███        ░███
!  █████     █████ ██████████░░█████████  █████   █████ █████ ░░░███████░      ░░░         ░░░
! ░░░░░     ░░░░░ ░░░░░░░░░░  ░░░░░░░░░  ░░░░░   ░░░░░ ░░░░░    ░░░░░░░
!
!
!  License:         MIT License
!                   meshio++ default license: LICENSE
!
!  Main authors:    Vicente Mataix Ferrandiz
!
!
! The meshio++ Fortran interface: a modern OO Fortran 2008 module layered on
! the C API (bindings_c/include/meshioplusplus/meshioplusplus.h) via
! ISO_C_BINDING, in the HDF5/PETSc style:
!
!     use meshioplusplus
!     type(mio_mesh) :: m
!     call m%read("bracket.msh")
!     print *, m%num_points()
!     call m%write("bracket.vtu")
!     call m%free()
!
! Conventions (differences from the C API):
!  - Arrays are Fortran-shaped: points are `points(dim, num_points)` and
!    connectivity `conn(nodes_per_cell, num_cells)`. Because Fortran is
!    column-major and the C core row-major, this is the SAME memory -- no
!    transpose happens anywhere in this module.
!  - Connectivity is 1-based here; the +-1 shift happens inside the copying
!    setters/getters (where a copy is made anyway). Zero-copy borrows
!    (points_ptr) carry no indices, so nothing zero-copy ever needs shifting.
!  - All indices (cell blocks, data names) are 1-based.
!  - Every fallible procedure takes `optional` `stat` (integer, 0 = success)
!    and `errmsg` (deferred-length character) arguments. If `stat` is absent
!    and the call fails, the message is printed and the program error stops
!    -- pass `stat` to handle errors yourself (the stdlib pattern).
!  - Handles are freed explicitly with `call m%free()` (no finalizer).
!
! Compiled .mod files are compiler-(major-version-)specific; this source file
! is installed next to the .mod so consumers on a different compiler can
! simply recompile the module (the HDF5 approach).
module meshioplusplus
    use, intrinsic :: iso_c_binding
    use, intrinsic :: iso_fortran_env, only: real32, real64, int32, int64, error_unit
    implicit none
    private

    public :: mio_mesh
    public :: mio_convert, mio_version, mio_mesh_backend, mio_error_message
    public :: mio_format_readable, mio_format_writable
    public :: mio_sniff_format

    ! mio_dtype values (must match the C enum in meshioplusplus.h).
    integer(c_int), parameter :: MIO_FLOAT32 = 0, MIO_FLOAT64 = 1
    integer(c_int), parameter :: MIO_INT32 = 4, MIO_INT64 = 5

    integer, parameter :: MIO_MAX_NDIM = 8
    integer, parameter :: STRBUF_LEN = 4096

    type :: mio_mesh
        private
        type(c_ptr) :: handle = c_null_ptr
    contains
        procedure :: create => mesh_create
        procedure :: free => mesh_free
        procedure :: is_valid => mesh_is_valid
        procedure :: read => mesh_read
        procedure :: write => mesh_write
        ! -- operations --
        procedure :: extract_surface => mesh_extract_surface
        procedure :: extract_skin => mesh_extract_skin
        procedure :: attach_quality => mesh_attach_quality
        procedure :: quality_counts => mesh_quality_counts
        ! -- building --
        procedure :: set_points => mesh_set_points
        procedure, private :: mesh_add_cell_block_i32
        procedure, private :: mesh_add_cell_block_i64
        generic :: add_cell_block => mesh_add_cell_block_i32, mesh_add_cell_block_i64
        procedure, private :: mesh_add_point_data_r1
        procedure, private :: mesh_add_point_data_r2
        generic :: add_point_data => mesh_add_point_data_r1, mesh_add_point_data_r2
        procedure :: add_cell_data => mesh_add_cell_data_r1
        procedure :: add_field_data => mesh_add_field_data_r1
        ! -- inspecting --
        procedure :: num_points => mesh_num_points
        procedure :: point_dim => mesh_point_dim
        procedure :: num_cell_blocks => mesh_num_cell_blocks
        procedure :: get_points => mesh_get_points
        procedure :: points_ptr => mesh_points_ptr
        procedure :: cell_block_type => mesh_cell_block_type
        procedure :: cell_block_num_cells => mesh_cell_block_num_cells
        procedure :: cell_block_nodes_per_cell => mesh_cell_block_nodes_per_cell
        procedure :: cell_block_is_ragged => mesh_cell_block_is_ragged
        procedure :: get_cell_block => mesh_get_cell_block
        procedure :: num_point_data => mesh_num_point_data
        procedure :: point_data_name => mesh_point_data_name
        procedure, private :: mesh_get_point_data_r1
        procedure, private :: mesh_get_point_data_r2
        generic :: get_point_data => mesh_get_point_data_r1, mesh_get_point_data_r2
        procedure :: num_cell_data => mesh_num_cell_data
        procedure :: cell_data_name => mesh_cell_data_name
        procedure :: cell_data_num_blocks => mesh_cell_data_num_blocks
        procedure :: get_cell_data => mesh_get_cell_data_r1
        procedure :: num_field_data => mesh_num_field_data
        procedure :: field_data_name => mesh_field_data_name
        procedure :: get_field_data => mesh_get_field_data_r1
    end type mio_mesh

    ! ------------------------------------------------------------------
    ! Raw bind(c) interfaces to libmeshioplusplus (private; the OO layer
    ! above is the public surface). Data pointers cross as type(c_ptr) via
    ! c_loc() so one C symbol serves every Fortran type/kind.
    ! ------------------------------------------------------------------
    interface
        function c_mio_version() bind(c, name="mio_version") result(p)
            import :: c_ptr
            type(c_ptr) :: p
        end function

        function c_mio_mesh_backend() bind(c, name="mio_mesh_backend") result(p)
            import :: c_ptr
            type(c_ptr) :: p
        end function

        function c_mio_last_error() bind(c, name="mio_last_error") result(p)
            import :: c_ptr
            type(c_ptr) :: p
        end function

        function c_mio_format_readable(format) bind(c, name="mio_format_readable") result(r)
            import :: c_char, c_int
            character(kind=c_char), dimension(*), intent(in) :: format
            integer(c_int) :: r
        end function

        function c_mio_format_writable(format) bind(c, name="mio_format_writable") result(r)
            import :: c_char, c_int
            character(kind=c_char), dimension(*), intent(in) :: format
            integer(c_int) :: r
        end function

        function c_mio_mesh_create() bind(c, name="mio_mesh_create") result(h)
            import :: c_ptr
            type(c_ptr) :: h
        end function

        subroutine c_mio_mesh_free(h) bind(c, name="mio_mesh_free")
            import :: c_ptr
            type(c_ptr), value :: h
        end subroutine

        function c_mio_read(path, format) bind(c, name="mio_read") result(h)
            import :: c_ptr, c_char
            character(kind=c_char), dimension(*), intent(in) :: path, format
            type(c_ptr) :: h
        end function

        function c_mio_write(path, h, format) bind(c, name="mio_write") result(s)
            import :: c_ptr, c_char, c_int
            character(kind=c_char), dimension(*), intent(in) :: path, format
            type(c_ptr), value :: h
            integer(c_int) :: s
        end function

        function c_mio_convert(in_path, in_format, out_path, out_format) &
                bind(c, name="mio_convert") result(s)
            import :: c_char, c_int
            character(kind=c_char), dimension(*), intent(in) :: in_path, in_format
            character(kind=c_char), dimension(*), intent(in) :: out_path, out_format
            integer(c_int) :: s
        end function

        function c_mio_extract_surface(h, record_parent_ids) &
                bind(c, name="mio_extract_surface") result(r)
            import :: c_ptr, c_int
            type(c_ptr), value :: h
            integer(c_int), value :: record_parent_ids
            type(c_ptr) :: r
        end function

        function c_mio_extract_skin(h, linearize) bind(c, name="mio_extract_skin") result(r)
            import :: c_ptr, c_int
            type(c_ptr), value :: h
            integer(c_int), value :: linearize
            type(c_ptr) :: r
        end function

        function c_mio_attach_quality(h) bind(c, name="mio_attach_quality") result(r)
            import :: c_ptr
            type(c_ptr), value :: h
            type(c_ptr) :: r
        end function

        function c_mio_quality_counts(h, nc, ninv, ndeg) &
                bind(c, name="mio_quality_counts") result(s)
            import :: c_ptr, c_int, c_int64_t
            type(c_ptr), value :: h
            integer(c_int64_t), intent(out) :: nc, ninv, ndeg
            integer(c_int) :: s
        end function

        function c_mio_sniff_format(path, buf, buflen) &
                bind(c, name="mio_sniff_format") result(n)
            import :: c_char, c_int64_t
            character(kind=c_char), dimension(*), intent(in) :: path
            character(kind=c_char), dimension(*), intent(inout) :: buf
            integer(c_int64_t), value :: buflen
            integer(c_int64_t) :: n
        end function

        function c_mio_mesh_set_points(h, dtype, num_points, dim, xyz) &
                bind(c, name="mio_mesh_set_points") result(s)
            import :: c_ptr, c_int, c_int64_t
            type(c_ptr), value :: h, xyz
            integer(c_int), value :: dtype
            integer(c_int64_t), value :: num_points, dim
            integer(c_int) :: s
        end function

        function c_mio_mesh_add_cell_block(h, cell_type, num_cells, nodes_per_cell, dtype, conn) &
                bind(c, name="mio_mesh_add_cell_block") result(s)
            import :: c_ptr, c_char, c_int, c_int64_t
            type(c_ptr), value :: h, conn
            character(kind=c_char), dimension(*), intent(in) :: cell_type
            integer(c_int64_t), value :: num_cells, nodes_per_cell
            integer(c_int), value :: dtype
            integer(c_int) :: s
        end function

        function c_mio_mesh_add_point_data(h, name, dtype, ndim, shape, data) &
                bind(c, name="mio_mesh_add_point_data") result(s)
            import :: c_ptr, c_char, c_int, c_int32_t, c_int64_t
            type(c_ptr), value :: h, data
            character(kind=c_char), dimension(*), intent(in) :: name
            integer(c_int), value :: dtype
            integer(c_int32_t), value :: ndim
            integer(c_int64_t), dimension(*), intent(in) :: shape
            integer(c_int) :: s
        end function

        function c_mio_mesh_append_cell_data(h, name, dtype, ndim, shape, data) &
                bind(c, name="mio_mesh_append_cell_data") result(s)
            import :: c_ptr, c_char, c_int, c_int32_t, c_int64_t
            type(c_ptr), value :: h, data
            character(kind=c_char), dimension(*), intent(in) :: name
            integer(c_int), value :: dtype
            integer(c_int32_t), value :: ndim
            integer(c_int64_t), dimension(*), intent(in) :: shape
            integer(c_int) :: s
        end function

        function c_mio_mesh_add_field_data(h, name, dtype, ndim, shape, data) &
                bind(c, name="mio_mesh_add_field_data") result(s)
            import :: c_ptr, c_char, c_int, c_int32_t, c_int64_t
            type(c_ptr), value :: h, data
            character(kind=c_char), dimension(*), intent(in) :: name
            integer(c_int), value :: dtype
            integer(c_int32_t), value :: ndim
            integer(c_int64_t), dimension(*), intent(in) :: shape
            integer(c_int) :: s
        end function

        function c_mio_mesh_num_points(h) bind(c, name="mio_mesh_num_points") result(n)
            import :: c_ptr, c_int64_t
            type(c_ptr), value :: h
            integer(c_int64_t) :: n
        end function

        function c_mio_mesh_point_dim(h) bind(c, name="mio_mesh_point_dim") result(n)
            import :: c_ptr, c_int64_t
            type(c_ptr), value :: h
            integer(c_int64_t) :: n
        end function

        function c_mio_mesh_get_points(h, data, dtype) &
                bind(c, name="mio_mesh_get_points") result(s)
            import :: c_ptr, c_int
            type(c_ptr), value :: h
            type(c_ptr), intent(out) :: data
            integer(c_int), intent(out) :: dtype
            integer(c_int) :: s
        end function

        function c_mio_mesh_num_cell_blocks(h) bind(c, name="mio_mesh_num_cell_blocks") result(n)
            import :: c_ptr, c_int64_t
            type(c_ptr), value :: h
            integer(c_int64_t) :: n
        end function

        function c_mio_mesh_cell_block_info(h, block, num_cells, nodes_per_cell, is_ragged) &
                bind(c, name="mio_mesh_cell_block_info") result(s)
            import :: c_ptr, c_int, c_int32_t, c_int64_t
            type(c_ptr), value :: h
            integer(c_int64_t), value :: block
            integer(c_int64_t), intent(out) :: num_cells, nodes_per_cell
            integer(c_int32_t), intent(out) :: is_ragged
            integer(c_int) :: s
        end function

        function c_mio_mesh_cell_block_type(h, block, buf, buflen) &
                bind(c, name="mio_mesh_cell_block_type") result(n)
            import :: c_ptr, c_char, c_int64_t
            type(c_ptr), value :: h
            integer(c_int64_t), value :: block, buflen
            character(kind=c_char), dimension(*), intent(inout) :: buf
            integer(c_int64_t) :: n
        end function

        function c_mio_mesh_cell_block_conn(h, block, conn, dtype) &
                bind(c, name="mio_mesh_cell_block_conn") result(s)
            import :: c_ptr, c_int, c_int64_t
            type(c_ptr), value :: h
            integer(c_int64_t), value :: block
            type(c_ptr), intent(out) :: conn
            integer(c_int), intent(out) :: dtype
            integer(c_int) :: s
        end function

        function c_mio_mesh_num_point_data(h) bind(c, name="mio_mesh_num_point_data") result(n)
            import :: c_ptr, c_int64_t
            type(c_ptr), value :: h
            integer(c_int64_t) :: n
        end function

        function c_mio_mesh_point_data_name(h, index, buf, buflen) &
                bind(c, name="mio_mesh_point_data_name") result(n)
            import :: c_ptr, c_char, c_int64_t
            type(c_ptr), value :: h
            integer(c_int64_t), value :: index, buflen
            character(kind=c_char), dimension(*), intent(inout) :: buf
            integer(c_int64_t) :: n
        end function

        function c_mio_mesh_get_point_data(h, name, data, dtype, ndim, shape) &
                bind(c, name="mio_mesh_get_point_data") result(s)
            import :: c_ptr, c_char, c_int, c_int32_t, c_int64_t
            type(c_ptr), value :: h
            character(kind=c_char), dimension(*), intent(in) :: name
            type(c_ptr), intent(out) :: data
            integer(c_int), intent(out) :: dtype
            integer(c_int32_t), intent(out) :: ndim
            integer(c_int64_t), dimension(*), intent(out) :: shape
            integer(c_int) :: s
        end function

        function c_mio_mesh_num_cell_data(h) bind(c, name="mio_mesh_num_cell_data") result(n)
            import :: c_ptr, c_int64_t
            type(c_ptr), value :: h
            integer(c_int64_t) :: n
        end function

        function c_mio_mesh_cell_data_name(h, index, buf, buflen) &
                bind(c, name="mio_mesh_cell_data_name") result(n)
            import :: c_ptr, c_char, c_int64_t
            type(c_ptr), value :: h
            integer(c_int64_t), value :: index, buflen
            character(kind=c_char), dimension(*), intent(inout) :: buf
            integer(c_int64_t) :: n
        end function

        function c_mio_mesh_cell_data_num_blocks(h, name) &
                bind(c, name="mio_mesh_cell_data_num_blocks") result(n)
            import :: c_ptr, c_char, c_int64_t
            type(c_ptr), value :: h
            character(kind=c_char), dimension(*), intent(in) :: name
            integer(c_int64_t) :: n
        end function

        function c_mio_mesh_get_cell_data(h, name, block, data, dtype, ndim, shape) &
                bind(c, name="mio_mesh_get_cell_data") result(s)
            import :: c_ptr, c_char, c_int, c_int32_t, c_int64_t
            type(c_ptr), value :: h
            character(kind=c_char), dimension(*), intent(in) :: name
            integer(c_int64_t), value :: block
            type(c_ptr), intent(out) :: data
            integer(c_int), intent(out) :: dtype
            integer(c_int32_t), intent(out) :: ndim
            integer(c_int64_t), dimension(*), intent(out) :: shape
            integer(c_int) :: s
        end function

        function c_mio_mesh_num_field_data(h) bind(c, name="mio_mesh_num_field_data") result(n)
            import :: c_ptr, c_int64_t
            type(c_ptr), value :: h
            integer(c_int64_t) :: n
        end function

        function c_mio_mesh_field_data_name(h, index, buf, buflen) &
                bind(c, name="mio_mesh_field_data_name") result(n)
            import :: c_ptr, c_char, c_int64_t
            type(c_ptr), value :: h
            integer(c_int64_t), value :: index, buflen
            character(kind=c_char), dimension(*), intent(inout) :: buf
            integer(c_int64_t) :: n
        end function

        function c_mio_mesh_get_field_data(h, name, data, dtype, ndim, shape) &
                bind(c, name="mio_mesh_get_field_data") result(s)
            import :: c_ptr, c_char, c_int, c_int32_t, c_int64_t
            type(c_ptr), value :: h
            character(kind=c_char), dimension(*), intent(in) :: name
            type(c_ptr), intent(out) :: data
            integer(c_int), intent(out) :: dtype
            integer(c_int32_t), intent(out) :: ndim
            integer(c_int64_t), dimension(*), intent(out) :: shape
            integer(c_int) :: s
        end function

        function c_strlen(s) bind(c, name="strlen") result(n)
            import :: c_ptr, c_size_t
            type(c_ptr), value :: s
            integer(c_size_t) :: n
        end function
    end interface

contains

    ! ------------------------------------------------------------------
    ! String / error helpers
    ! ------------------------------------------------------------------

    !> NUL-terminated copy for passing to C ("" stays "" + NUL, meaning
    !> "infer format" on the C side).
    pure function c_str(f) result(c)
        character(*), intent(in) :: f
        character(kind=c_char, len=:), allocatable :: c
        c = trim(f)//c_null_char
    end function

    !> Fortran string from a NUL-terminated C pointer (static storage).
    function c_ptr_to_string(p) result(s)
        type(c_ptr), intent(in) :: p
        character(:), allocatable :: s
        character(kind=c_char), pointer :: chars(:)
        integer :: n, i
        s = ''
        if (.not. c_associated(p)) return
        n = int(c_strlen(p))
        if (n <= 0) return
        call c_f_pointer(p, chars, [n])
        s = repeat(' ', n)  ! reallocation on assignment (s starts as '')
        do i = 1, n
            s(i:i) = chars(i)
        end do
    end function

    !> Fortran string from the first `n` chars of a C char buffer.
    function from_c_buf(buf, n) result(s)
        character(kind=c_char), intent(in) :: buf(*)
        integer, intent(in) :: n
        character(:), allocatable :: s
        integer :: i
        allocate (character(max(n, 0)) :: s)
        do i = 1, n
            s(i:i) = buf(i)
        end do
    end function

    !> The failure message of the most recent failed meshio++ call on this
    !> thread ('' if none).
    function mio_error_message() result(msg)
        character(:), allocatable :: msg
        msg = c_ptr_to_string(c_mio_last_error())
    end function

    !> Map a C status to the optional stat/errmsg pair; no stat + failure =
    !> print and error stop.
    subroutine handle_status(status, what, stat, errmsg)
        integer(c_int), intent(in) :: status
        character(*), intent(in) :: what
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        if (present(stat)) stat = int(status)
        if (present(errmsg)) errmsg = ''
        if (status /= 0_c_int) then
            if (present(errmsg)) errmsg = mio_error_message()
            if (.not. present(stat)) then
                write (error_unit, '(a)') 'meshio++ ('//what//'): '//mio_error_message()
                error stop 1
            end if
        end if
    end subroutine

    !> Report a Fortran-side failure (bad handle, shape mismatch, ...)
    !> through the same stat/errmsg protocol.
    subroutine handle_failure(what, msg, stat, errmsg)
        character(*), intent(in) :: what, msg
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        if (present(stat)) stat = 1
        if (present(errmsg)) errmsg = msg
        if (.not. present(stat)) then
            write (error_unit, '(a)') 'meshio++ ('//what//'): '//msg
            error stop 1
        end if
    end subroutine

    subroutine clear_status(stat, errmsg)
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        if (present(stat)) stat = 0
        if (present(errmsg)) errmsg = ''
    end subroutine

    ! ------------------------------------------------------------------
    ! Module-level procedures
    ! ------------------------------------------------------------------

    !> The meshio++ version string, e.g. "6.1.0".
    function mio_version() result(v)
        character(:), allocatable :: v
        v = c_ptr_to_string(c_mio_version())
    end function

    !> The compile-time mesh backend: "meshio", "native", or "kratos".
    function mio_mesh_backend() result(b)
        character(:), allocatable :: b
        b = c_ptr_to_string(c_mio_mesh_backend())
    end function

    !> .true. if `format` (e.g. "gmsh", "vtu", "med") is readable in this build.
    function mio_format_readable(format) result(r)
        character(*), intent(in) :: format
        logical :: r
        r = c_mio_format_readable(c_str(format)) /= 0_c_int
    end function

    !> .true. if `format` is writable in this build.
    function mio_format_writable(format) result(r)
        character(*), intent(in) :: format
        logical :: r
        r = c_mio_format_writable(c_str(format)) /= 0_c_int
    end function

    !> Read `in_path` and immediately write it to `out_path` (the CLI's
    !> `convert`). Formats are inferred from the extensions unless given.
    subroutine mio_convert(in_path, out_path, in_format, out_format, stat, errmsg)
        character(*), intent(in) :: in_path, out_path
        character(*), intent(in), optional :: in_format, out_format
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        character(:), allocatable :: ifmt, ofmt
        ifmt = ''; if (present(in_format)) ifmt = in_format
        ofmt = ''; if (present(out_format)) ofmt = out_format
        call handle_status(c_mio_convert(c_str(in_path), c_str(ifmt), c_str(out_path), &
                                         c_str(ofmt)), 'convert', stat, errmsg)
    end subroutine

    !> Guess a mesh file's format from its contents ('' if undetermined).
    function mio_sniff_format(path) result(fmt)
        character(*), intent(in) :: path
        character(:), allocatable :: fmt
        character(kind=c_char) :: buf(STRBUF_LEN)
        integer(c_int64_t) :: n
        fmt = ''
        n = c_mio_sniff_format(c_str(path), buf, int(STRBUF_LEN, c_int64_t))
        if (n > 0) fmt = from_c_buf(buf, min(int(n), STRBUF_LEN - 1))
    end function

    ! ------------------------------------------------------------------
    ! mio_mesh: lifecycle & file I/O
    ! ------------------------------------------------------------------

    !> Allocate an empty mesh (read() does this implicitly).
    subroutine mesh_create(self, stat, errmsg)
        class(mio_mesh), intent(inout) :: self
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        call mesh_free(self)
        self%handle = c_mio_mesh_create()
        if (.not. c_associated(self%handle)) then
            call handle_failure('create', mio_error_message(), stat, errmsg)
            return
        end if
        call clear_status(stat, errmsg)
    end subroutine

    !> Release the mesh and every pointer borrowed from it. Idempotent.
    subroutine mesh_free(self)
        class(mio_mesh), intent(inout) :: self
        if (c_associated(self%handle)) call c_mio_mesh_free(self%handle)
        self%handle = c_null_ptr
    end subroutine

    !> .true. between a successful create()/read() and free().
    logical function mesh_is_valid(self)
        class(mio_mesh), intent(in) :: self
        mesh_is_valid = c_associated(self%handle)
    end function

    !> Read a mesh file, replacing any previous content of this handle.
    subroutine mesh_read(self, path, format, stat, errmsg)
        class(mio_mesh), intent(inout) :: self
        character(*), intent(in) :: path
        character(*), intent(in), optional :: format
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        character(:), allocatable :: fmt
        type(c_ptr) :: h
        fmt = ''; if (present(format)) fmt = format
        h = c_mio_read(c_str(path), c_str(fmt))
        if (.not. c_associated(h)) then
            call handle_failure('read', mio_error_message(), stat, errmsg)
            return
        end if
        call mesh_free(self)
        self%handle = h
        call clear_status(stat, errmsg)
    end subroutine

    !> Write the mesh to a file.
    subroutine mesh_write(self, path, format, stat, errmsg)
        class(mio_mesh), intent(in) :: self
        character(*), intent(in) :: path
        character(*), intent(in), optional :: format
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        character(:), allocatable :: fmt
        fmt = ''; if (present(format)) fmt = format
        call handle_status(c_mio_write(c_str(path), self%handle, c_str(fmt)), 'write', &
                           stat, errmsg)
    end subroutine

    !> Extract the boundary of the mesh's highest-dimension cells as a new mesh
    !> (volume -> faces, surface -> edges).
    function mesh_extract_surface(self, record_parent_ids, stat, errmsg) result(out)
        class(mio_mesh), intent(in) :: self
        logical, intent(in), optional :: record_parent_ids
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(mio_mesh) :: out
        integer(c_int) :: rec
        rec = 0
        if (present(record_parent_ids)) then
            if (record_parent_ids) rec = 1
        end if
        out%handle = c_mio_extract_surface(self%handle, rec)
        if (.not. c_associated(out%handle)) then
            call handle_failure('extract_surface', mio_error_message(), stat, errmsg)
            return
        end if
        call clear_status(stat, errmsg)
    end function

    !> Extract the boundary skin of a volume mesh as a new surface mesh.
    function mesh_extract_skin(self, linearize, stat, errmsg) result(out)
        class(mio_mesh), intent(in) :: self
        logical, intent(in), optional :: linearize
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(mio_mesh) :: out
        integer(c_int) :: lin
        lin = 0
        if (present(linearize)) then
            if (linearize) lin = 1
        end if
        out%handle = c_mio_extract_skin(self%handle, lin)
        if (.not. c_associated(out%handle)) then
            call handle_failure('extract_skin', mio_error_message(), stat, errmsg)
            return
        end if
        call clear_status(stat, errmsg)
    end function

    !> A copy of the mesh with per-cell quality metrics attached as cell_data.
    function mesh_attach_quality(self, stat, errmsg) result(out)
        class(mio_mesh), intent(in) :: self
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(mio_mesh) :: out
        out%handle = c_mio_attach_quality(self%handle)
        if (.not. c_associated(out%handle)) then
            call handle_failure('attach_quality', mio_error_message(), stat, errmsg)
            return
        end if
        call clear_status(stat, errmsg)
    end function

    !> Aggregate quality counts (total / inverted / degenerate cells).
    subroutine mesh_quality_counts(self, num_cells, num_inverted, num_degenerate, stat, errmsg)
        class(mio_mesh), intent(in) :: self
        integer(int64), intent(out), optional :: num_cells, num_inverted, num_degenerate
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        integer(c_int64_t) :: nc, ninv, ndeg
        if (c_mio_quality_counts(self%handle, nc, ninv, ndeg) /= 0_c_int) then
            call handle_failure('quality_counts', mio_error_message(), stat, errmsg)
            return
        end if
        if (present(num_cells)) num_cells = nc
        if (present(num_inverted)) num_inverted = ninv
        if (present(num_degenerate)) num_degenerate = ndeg
        call clear_status(stat, errmsg)
    end subroutine

    ! ------------------------------------------------------------------
    ! mio_mesh: building (setters copy; see module header for layout rules)
    ! ------------------------------------------------------------------

    !> Assign the point coordinates from a `points(dim, num_points)` array.
    subroutine mesh_set_points(self, points, stat, errmsg)
        class(mio_mesh), intent(inout) :: self
        real(real64), intent(in), contiguous, target :: points(:, :)
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        call ensure_handle(self, stat, errmsg)
        if (.not. c_associated(self%handle)) return
        call handle_status(c_mio_mesh_set_points(self%handle, MIO_FLOAT64, &
                                                 int(size(points, 2), c_int64_t), &
                                                 int(size(points, 1), c_int64_t), &
                                                 c_loc(points)), 'set_points', stat, errmsg)
    end subroutine

    !> Append one cell block from a 1-based `conn(nodes_per_cell, num_cells)`
    !> array (default-integer version).
    subroutine mesh_add_cell_block_i32(self, cell_type, conn, stat, errmsg)
        class(mio_mesh), intent(inout) :: self
        character(*), intent(in) :: cell_type
        integer(int32), intent(in) :: conn(:, :)
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        integer(c_int64_t), allocatable, target :: shifted(:, :)
        call ensure_handle(self, stat, errmsg)
        if (.not. c_associated(self%handle)) return
        shifted = int(conn, c_int64_t) - 1_c_int64_t
        call handle_status(c_mio_mesh_add_cell_block(self%handle, c_str(cell_type), &
                                                     int(size(conn, 2), c_int64_t), &
                                                     int(size(conn, 1), c_int64_t), MIO_INT64, &
                                                     c_loc(shifted)), 'add_cell_block', &
                           stat, errmsg)
    end subroutine

    !> As above, for `integer(int64)` connectivity.
    subroutine mesh_add_cell_block_i64(self, cell_type, conn, stat, errmsg)
        class(mio_mesh), intent(inout) :: self
        character(*), intent(in) :: cell_type
        integer(int64), intent(in) :: conn(:, :)
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        integer(c_int64_t), allocatable, target :: shifted(:, :)
        call ensure_handle(self, stat, errmsg)
        if (.not. c_associated(self%handle)) return
        shifted = int(conn, c_int64_t) - 1_c_int64_t
        call handle_status(c_mio_mesh_add_cell_block(self%handle, c_str(cell_type), &
                                                     int(size(conn, 2), c_int64_t), &
                                                     int(size(conn, 1), c_int64_t), MIO_INT64, &
                                                     c_loc(shifted)), 'add_cell_block', &
                           stat, errmsg)
    end subroutine

    !> Attach a scalar per-point field: `data(num_points)`.
    subroutine mesh_add_point_data_r1(self, name, data, stat, errmsg)
        class(mio_mesh), intent(inout) :: self
        character(*), intent(in) :: name
        real(real64), intent(in), contiguous, target :: data(:)
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        integer(c_int64_t) :: shape(1)
        call ensure_handle(self, stat, errmsg)
        if (.not. c_associated(self%handle)) return
        shape(1) = int(size(data), c_int64_t)
        call handle_status(c_mio_mesh_add_point_data(self%handle, c_str(name), MIO_FLOAT64, &
                                                     1_c_int32_t, shape, c_loc(data)), &
                           'add_point_data', stat, errmsg)
    end subroutine

    !> Attach a vector per-point field: `data(num_components, num_points)`
    !> (same memory as the C API's row-major `(num_points, num_components)`).
    subroutine mesh_add_point_data_r2(self, name, data, stat, errmsg)
        class(mio_mesh), intent(inout) :: self
        character(*), intent(in) :: name
        real(real64), intent(in), contiguous, target :: data(:, :)
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        integer(c_int64_t) :: shape(2)
        call ensure_handle(self, stat, errmsg)
        if (.not. c_associated(self%handle)) return
        shape(1) = int(size(data, 2), c_int64_t)  ! C shape is the reverse of the Fortran one
        shape(2) = int(size(data, 1), c_int64_t)
        call handle_status(c_mio_mesh_add_point_data(self%handle, c_str(name), MIO_FLOAT64, &
                                                     2_c_int32_t, shape, c_loc(data)), &
                           'add_point_data', stat, errmsg)
    end subroutine

    !> Append the named cell-data field's array for the next cell block (call
    !> once per block, in block order): `data(num_cells_in_block)`.
    subroutine mesh_add_cell_data_r1(self, name, data, stat, errmsg)
        class(mio_mesh), intent(inout) :: self
        character(*), intent(in) :: name
        real(real64), intent(in), contiguous, target :: data(:)
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        integer(c_int64_t) :: shape(1)
        call ensure_handle(self, stat, errmsg)
        if (.not. c_associated(self%handle)) return
        shape(1) = int(size(data), c_int64_t)
        call handle_status(c_mio_mesh_append_cell_data(self%handle, c_str(name), MIO_FLOAT64, &
                                                       1_c_int32_t, shape, c_loc(data)), &
                           'add_cell_data', stat, errmsg)
    end subroutine

    !> Attach a named mesh-level (field) data array.
    subroutine mesh_add_field_data_r1(self, name, data, stat, errmsg)
        class(mio_mesh), intent(inout) :: self
        character(*), intent(in) :: name
        real(real64), intent(in), contiguous, target :: data(:)
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        integer(c_int64_t) :: shape(1)
        call ensure_handle(self, stat, errmsg)
        if (.not. c_associated(self%handle)) return
        shape(1) = int(size(data), c_int64_t)
        call handle_status(c_mio_mesh_add_field_data(self%handle, c_str(name), MIO_FLOAT64, &
                                                     1_c_int32_t, shape, c_loc(data)), &
                           'add_field_data', stat, errmsg)
    end subroutine

    ! ------------------------------------------------------------------
    ! mio_mesh: inspection
    ! ------------------------------------------------------------------

    integer(int64) function mesh_num_points(self)
        class(mio_mesh), intent(in) :: self
        mesh_num_points = c_mio_mesh_num_points(self%handle)
    end function

    integer(int64) function mesh_point_dim(self)
        class(mio_mesh), intent(in) :: self
        mesh_point_dim = c_mio_mesh_point_dim(self%handle)
    end function

    integer(int64) function mesh_num_cell_blocks(self)
        class(mio_mesh), intent(in) :: self
        mesh_num_cell_blocks = c_mio_mesh_num_cell_blocks(self%handle)
    end function

    !> Copy the coordinates into `points(point_dim, num_points)` (allocated
    !> here), converting to real64 if the mesh stores single precision.
    subroutine mesh_get_points(self, points, stat, errmsg)
        class(mio_mesh), intent(in) :: self
        real(real64), allocatable, intent(out) :: points(:, :)
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(c_ptr) :: p
        integer(c_int) :: dtype, status
        integer(int64) :: n, dim
        real(c_double), pointer :: d64(:)
        real(c_float), pointer :: d32(:)
        status = c_mio_mesh_get_points(self%handle, p, dtype)
        if (status /= 0_c_int) then
            call handle_status(status, 'get_points', stat, errmsg)
            return
        end if
        n = c_mio_mesh_num_points(self%handle)
        dim = c_mio_mesh_point_dim(self%handle)
        allocate (points(dim, n))
        if (n*dim == 0) then
            call clear_status(stat, errmsg)
            return
        end if
        select case (dtype)
        case (MIO_FLOAT64)
            call c_f_pointer(p, d64, [n*dim])
            points = reshape(d64, [dim, n])
        case (MIO_FLOAT32)
            call c_f_pointer(p, d32, [n*dim])
            points = reshape(real(d32, real64), [dim, n])
        case default
            call handle_failure('get_points', 'unexpected points dtype', stat, errmsg)
            return
        end select
        call clear_status(stat, errmsg)
    end subroutine

    !> Zero-copy borrow of the coordinates as `ptr(point_dim, num_points)`.
    !> Valid until the next mutating call on this mesh or free(); fails if
    !> the mesh stores anything but real64 (then use get_points).
    subroutine mesh_points_ptr(self, ptr, stat, errmsg)
        class(mio_mesh), intent(in) :: self
        real(real64), pointer, intent(out) :: ptr(:, :)
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(c_ptr) :: p
        integer(c_int) :: dtype, status
        ptr => null()
        status = c_mio_mesh_get_points(self%handle, p, dtype)
        if (status /= 0_c_int) then
            call handle_status(status, 'points_ptr', stat, errmsg)
            return
        end if
        if (dtype /= MIO_FLOAT64) then
            call handle_failure('points_ptr', 'points are not real64; use get_points', &
                                stat, errmsg)
            return
        end if
        call c_f_pointer(p, ptr, [c_mio_mesh_point_dim(self%handle), &
                                  c_mio_mesh_num_points(self%handle)])
        call clear_status(stat, errmsg)
    end subroutine

    !> The meshio type name of 1-based cell block `block` (e.g. "tetra10").
    function mesh_cell_block_type(self, block, stat, errmsg) result(type_name)
        class(mio_mesh), intent(in) :: self
        integer, intent(in) :: block
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        character(:), allocatable :: type_name
        character(kind=c_char) :: buf(STRBUF_LEN)
        integer(c_int64_t) :: n
        type_name = ''
        n = c_mio_mesh_cell_block_type(self%handle, int(block - 1, c_int64_t), buf, &
                                       int(STRBUF_LEN, c_int64_t))
        if (n < 0) then
            call handle_failure('cell_block_type', mio_error_message(), stat, errmsg)
            return
        end if
        type_name = from_c_buf(buf, min(int(n), STRBUF_LEN - 1))
        call clear_status(stat, errmsg)
    end function

    integer(int64) function mesh_cell_block_num_cells(self, block)
        class(mio_mesh), intent(in) :: self
        integer, intent(in) :: block
        integer(c_int64_t) :: nc, npc
        integer(c_int32_t) :: ragged
        mesh_cell_block_num_cells = -1
        if (c_mio_mesh_cell_block_info(self%handle, int(block - 1, c_int64_t), nc, npc, &
                                       ragged) == 0_c_int) mesh_cell_block_num_cells = nc
    end function

    integer(int64) function mesh_cell_block_nodes_per_cell(self, block)
        class(mio_mesh), intent(in) :: self
        integer, intent(in) :: block
        integer(c_int64_t) :: nc, npc
        integer(c_int32_t) :: ragged
        mesh_cell_block_nodes_per_cell = -1
        if (c_mio_mesh_cell_block_info(self%handle, int(block - 1, c_int64_t), nc, npc, &
                                       ragged) == 0_c_int) mesh_cell_block_nodes_per_cell = npc
    end function

    logical function mesh_cell_block_is_ragged(self, block)
        class(mio_mesh), intent(in) :: self
        integer, intent(in) :: block
        integer(c_int64_t) :: nc, npc
        integer(c_int32_t) :: ragged
        mesh_cell_block_is_ragged = .false.
        if (c_mio_mesh_cell_block_info(self%handle, int(block - 1, c_int64_t), nc, npc, &
                                       ragged) == 0_c_int) mesh_cell_block_is_ragged = ragged /= 0
    end function

    !> Copy 1-based cell block `block` into `conn(nodes_per_cell, num_cells)`
    !> (allocated here), shifting the 0-based core indices to 1-based.
    subroutine mesh_get_cell_block(self, block, conn, stat, errmsg)
        class(mio_mesh), intent(in) :: self
        integer, intent(in) :: block
        integer(int64), allocatable, intent(out) :: conn(:, :)
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(c_ptr) :: p
        integer(c_int) :: dtype, status
        integer(c_int64_t) :: nc, npc
        integer(c_int32_t) :: ragged
        integer(c_int64_t), pointer :: i64(:)
        integer(c_int32_t), pointer :: i32(:)
        status = c_mio_mesh_cell_block_info(self%handle, int(block - 1, c_int64_t), nc, npc, &
                                            ragged)
        if (status == 0_c_int) status = c_mio_mesh_cell_block_conn(self%handle, &
                                                                   int(block - 1, c_int64_t), &
                                                                   p, dtype)
        if (status /= 0_c_int) then
            call handle_status(status, 'get_cell_block', stat, errmsg)
            return
        end if
        allocate (conn(npc, nc))
        if (nc*npc == 0) then
            call clear_status(stat, errmsg)
            return
        end if
        select case (dtype)
        case (MIO_INT64)
            call c_f_pointer(p, i64, [nc*npc])
            conn = reshape(i64, [npc, nc]) + 1_int64
        case (MIO_INT32)
            call c_f_pointer(p, i32, [nc*npc])
            conn = reshape(int(i32, int64), [npc, nc]) + 1_int64
        case default
            call handle_failure('get_cell_block', 'unexpected connectivity dtype', stat, errmsg)
            return
        end select
        call clear_status(stat, errmsg)
    end subroutine

    integer(int64) function mesh_num_point_data(self)
        class(mio_mesh), intent(in) :: self
        mesh_num_point_data = c_mio_mesh_num_point_data(self%handle)
    end function

    !> The 1-based `index`-th point-data name (ascending lexicographic order).
    function mesh_point_data_name(self, index, stat, errmsg) result(name)
        class(mio_mesh), intent(in) :: self
        integer, intent(in) :: index
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        character(:), allocatable :: name
        character(kind=c_char) :: buf(STRBUF_LEN)
        integer(c_int64_t) :: n
        name = ''
        n = c_mio_mesh_point_data_name(self%handle, int(index - 1, c_int64_t), buf, &
                                       int(STRBUF_LEN, c_int64_t))
        if (n < 0) then
            call handle_failure('point_data_name', mio_error_message(), stat, errmsg)
            return
        end if
        name = from_c_buf(buf, min(int(n), STRBUF_LEN - 1))
        call clear_status(stat, errmsg)
    end function

    !> Copy the named scalar point-data array into `data(:)` (allocated here).
    subroutine mesh_get_point_data_r1(self, name, data, stat, errmsg)
        class(mio_mesh), intent(in) :: self
        character(*), intent(in) :: name
        real(real64), allocatable, intent(out) :: data(:)
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        call get_named_r1(self, 'point_data', name, -1_int64, data, stat, errmsg)
    end subroutine

    !> Copy the named vector point-data array into
    !> `data(num_components, num_points)` (allocated here).
    subroutine mesh_get_point_data_r2(self, name, data, stat, errmsg)
        class(mio_mesh), intent(in) :: self
        character(*), intent(in) :: name
        real(real64), allocatable, intent(out) :: data(:, :)
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(c_ptr) :: p
        integer(c_int) :: dtype, status
        integer(c_int32_t) :: ndim
        integer(c_int64_t) :: shape(MIO_MAX_NDIM)
        status = c_mio_mesh_get_point_data(self%handle, c_str(name), p, dtype, ndim, shape)
        if (status /= 0_c_int) then
            call handle_status(status, 'get_point_data', stat, errmsg)
            return
        end if
        if (ndim /= 2_c_int32_t) then
            call handle_failure('get_point_data', 'point_data "'//trim(name)// &
                                '" is not rank-2', stat, errmsg)
            return
        end if
        call copy_out_r2(p, dtype, shape(1), shape(2), data, 'get_point_data', stat, errmsg)
    end subroutine

    integer(int64) function mesh_num_cell_data(self)
        class(mio_mesh), intent(in) :: self
        mesh_num_cell_data = c_mio_mesh_num_cell_data(self%handle)
    end function

    !> The 1-based `index`-th cell-data name (sorted, as above).
    function mesh_cell_data_name(self, index, stat, errmsg) result(name)
        class(mio_mesh), intent(in) :: self
        integer, intent(in) :: index
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        character(:), allocatable :: name
        character(kind=c_char) :: buf(STRBUF_LEN)
        integer(c_int64_t) :: n
        name = ''
        n = c_mio_mesh_cell_data_name(self%handle, int(index - 1, c_int64_t), buf, &
                                      int(STRBUF_LEN, c_int64_t))
        if (n < 0) then
            call handle_failure('cell_data_name', mio_error_message(), stat, errmsg)
            return
        end if
        name = from_c_buf(buf, min(int(n), STRBUF_LEN - 1))
        call clear_status(stat, errmsg)
    end function

    integer(int64) function mesh_cell_data_num_blocks(self, name)
        class(mio_mesh), intent(in) :: self
        character(*), intent(in) :: name
        mesh_cell_data_num_blocks = c_mio_mesh_cell_data_num_blocks(self%handle, c_str(name))
    end function

    !> Copy the named cell-data field's array for 1-based cell block `block`
    !> into `data(:)` (allocated here).
    subroutine mesh_get_cell_data_r1(self, name, block, data, stat, errmsg)
        class(mio_mesh), intent(in) :: self
        character(*), intent(in) :: name
        integer, intent(in) :: block
        real(real64), allocatable, intent(out) :: data(:)
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        call get_named_r1(self, 'cell_data', name, int(block - 1, int64), data, stat, errmsg)
    end subroutine

    integer(int64) function mesh_num_field_data(self)
        class(mio_mesh), intent(in) :: self
        mesh_num_field_data = c_mio_mesh_num_field_data(self%handle)
    end function

    !> The 1-based `index`-th field-data name (sorted, as above).
    function mesh_field_data_name(self, index, stat, errmsg) result(name)
        class(mio_mesh), intent(in) :: self
        integer, intent(in) :: index
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        character(:), allocatable :: name
        character(kind=c_char) :: buf(STRBUF_LEN)
        integer(c_int64_t) :: n
        name = ''
        n = c_mio_mesh_field_data_name(self%handle, int(index - 1, c_int64_t), buf, &
                                       int(STRBUF_LEN, c_int64_t))
        if (n < 0) then
            call handle_failure('field_data_name', mio_error_message(), stat, errmsg)
            return
        end if
        name = from_c_buf(buf, min(int(n), STRBUF_LEN - 1))
        call clear_status(stat, errmsg)
    end function

    !> Copy the named field-data array into `data(:)` (allocated here).
    subroutine mesh_get_field_data_r1(self, name, data, stat, errmsg)
        class(mio_mesh), intent(in) :: self
        character(*), intent(in) :: name
        real(real64), allocatable, intent(out) :: data(:)
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        call get_named_r1(self, 'field_data', name, -2_int64, data, stat, errmsg)
    end subroutine

    ! ------------------------------------------------------------------
    ! Private implementation helpers
    ! ------------------------------------------------------------------

    !> Setters on a never-created handle allocate one on the fly, so
    !> `type(mio_mesh) :: m; call m%set_points(...)` just works.
    subroutine ensure_handle(self, stat, errmsg)
        class(mio_mesh), intent(inout) :: self
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        if (.not. c_associated(self%handle)) call mesh_create(self, stat, errmsg)
    end subroutine

    !> Shared rank-1 copy-getter over the three named-data families
    !> (block >= 0: cell_data block; -1: point_data; -2: field_data).
    subroutine get_named_r1(self, family, name, block, data, stat, errmsg)
        class(mio_mesh), intent(in) :: self
        character(*), intent(in) :: family, name
        integer(int64), intent(in) :: block
        real(real64), allocatable, intent(out) :: data(:)
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        type(c_ptr) :: p
        integer(c_int) :: dtype, status
        integer(c_int32_t) :: ndim
        integer(c_int64_t) :: shape(MIO_MAX_NDIM)
        select case (family)
        case ('point_data')
            status = c_mio_mesh_get_point_data(self%handle, c_str(name), p, dtype, ndim, shape)
        case ('cell_data')
            status = c_mio_mesh_get_cell_data(self%handle, c_str(name), block, p, dtype, ndim, &
                                              shape)
        case default
            status = c_mio_mesh_get_field_data(self%handle, c_str(name), p, dtype, ndim, shape)
        end select
        if (status /= 0_c_int) then
            call handle_status(status, 'get_'//family, stat, errmsg)
            return
        end if
        if (ndim /= 1_c_int32_t) then
            call handle_failure('get_'//family, family//' "'//trim(name)// &
                                '" is not rank-1 (use the rank-2 getter)', stat, errmsg)
            return
        end if
        call copy_out_r1(p, dtype, shape(1), data, 'get_'//family, stat, errmsg)
    end subroutine

    !> Copy a C array of `n` elements of any supported dtype into a real64
    !> rank-1 allocatable.
    subroutine copy_out_r1(p, dtype, n, data, what, stat, errmsg)
        type(c_ptr), intent(in) :: p
        integer(c_int), intent(in) :: dtype
        integer(c_int64_t), intent(in) :: n
        real(real64), allocatable, intent(out) :: data(:)
        character(*), intent(in) :: what
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        real(c_double), pointer :: d64(:)
        real(c_float), pointer :: d32(:)
        integer(c_int64_t), pointer :: i64(:)
        integer(c_int32_t), pointer :: i32(:)
        allocate (data(n))
        if (n == 0) then
            call clear_status(stat, errmsg)
            return
        end if
        select case (dtype)
        case (MIO_FLOAT64)
            call c_f_pointer(p, d64, [n]); data = d64
        case (MIO_FLOAT32)
            call c_f_pointer(p, d32, [n]); data = real(d32, real64)
        case (MIO_INT64)
            call c_f_pointer(p, i64, [n]); data = real(i64, real64)
        case (MIO_INT32)
            call c_f_pointer(p, i32, [n]); data = real(i32, real64)
        case default
            call handle_failure(what, 'unsupported dtype for real64 copy', stat, errmsg)
            return
        end select
        call clear_status(stat, errmsg)
    end subroutine

    !> As copy_out_r1 for a C row-major `(n0, n1)` array, delivered as the
    !> Fortran `(n1, n0)` view of the same memory order.
    subroutine copy_out_r2(p, dtype, n0, n1, data, what, stat, errmsg)
        type(c_ptr), intent(in) :: p
        integer(c_int), intent(in) :: dtype
        integer(c_int64_t), intent(in) :: n0, n1
        real(real64), allocatable, intent(out) :: data(:, :)
        character(*), intent(in) :: what
        integer, intent(out), optional :: stat
        character(:), allocatable, intent(out), optional :: errmsg
        real(c_double), pointer :: d64(:)
        real(c_float), pointer :: d32(:)
        integer(c_int64_t), pointer :: i64(:)
        integer(c_int32_t), pointer :: i32(:)
        allocate (data(n1, n0))
        if (n0*n1 == 0) then
            call clear_status(stat, errmsg)
            return
        end if
        select case (dtype)
        case (MIO_FLOAT64)
            call c_f_pointer(p, d64, [n0*n1]); data = reshape(d64, [n1, n0])
        case (MIO_FLOAT32)
            call c_f_pointer(p, d32, [n0*n1]); data = reshape(real(d32, real64), [n1, n0])
        case (MIO_INT64)
            call c_f_pointer(p, i64, [n0*n1]); data = reshape(real(i64, real64), [n1, n0])
        case (MIO_INT32)
            call c_f_pointer(p, i32, [n0*n1]); data = reshape(real(i32, real64), [n1, n0])
        case default
            call handle_failure(what, 'unsupported dtype for real64 copy', stat, errmsg)
            return
        end select
        call clear_status(stat, errmsg)
    end subroutine

end module meshioplusplus
