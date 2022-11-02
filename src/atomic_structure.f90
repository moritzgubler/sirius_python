module atomic_structure
    implicit none
    type geometry
        integer :: nat
        real(8), dimension(:, :), allocatable :: pos
        real(8), dimension(3, 3) :: lat
        character(len=2), dimension(:), allocatable :: atomnames
        contains

        procedure :: initialize
        procedure :: dealloc_arrays
        procedure :: broadcast
        procedure :: recv

    end type geometry
contains

    subroutine initialize(t, nat, pos, lat, atomnames)
        class(geometry) :: t
        integer :: nat
        real(8) :: pos(3, nat)
        real(8) :: lat(3,3)
        character(len=2), dimension(nat) :: atomnames

        if ( allocated(t%pos)) deallocate(t%pos)
        if ( allocated(t%atomnames) ) deallocate(t%atomnames)
        allocate(t%pos(3, nat))
        allocate(t%atomnames(nat))
        t%nat = nat
        t%pos = pos
        t%lat = lat
        t%atomnames = atomnames
    end subroutine initialize

    subroutine dealloc_arrays(t)
        class(geometry) :: t
        if ( allocated(t%pos) ) deallocate(t%pos)
        if ( allocated(t%atomnames) ) deallocate(t%atomnames)
    end subroutine dealloc_arrays

    subroutine broadcast(t)
        use mpi
        class(geometry) :: t
        integer :: ierr, rank
        character(len=2*t%nat) :: atoms_string


        call mpi_comm_rank(MPI_COMM_WORLD, rank, ierr)

        call string_array_to_string(atoms_string, t%atomnames, t%nat, 2)

        if ( .not. allocated(t%pos)) stop "position in geometry broadcast not initialized"
        if ( .not. allocated(t%atomnames)) stop 'atomnames in geometry broadcast not initialized'

        call mpi_bcast(t%nat, 1, MPI_INT, rank, MPI_COMM_WORLD, ierr)
        if (ierr /= 0) print*, 'sbjoi'
        call mpi_bcast(t%pos, 3*t%nat, MPI_DOUBLE, rank, MPI_COMM_WORLD, ierr)
        if (ierr /= 0) print*, 'sbjasdfg;o'
        call mpi_bcast(t%nat, 9, MPI_DOUBLE, rank, MPI_COMM_WORLD, ierr)
        if (ierr /= 0) print*, 'stionkwreg'
        call mpi_bcast(atoms_string, 2*t%nat, MPI_CHARACTER, rank, MPI_COMM_WORLD, ierr)
        if (ierr /= 0) print*, '67yu90os'

    contains
        subroutine string_array_to_string(str, str_arr, len_array, len_elem)
          implicit none
          integer, intent(in) :: len_array
          integer, intent(in) :: len_elem
          character(len=len_array*len_elem), intent(out) :: str
          character(len=len_elem), dimension(len_array), intent(in) :: str_arr
          integer :: i
          do i = 1, len_array, 1
            str((i*len_elem - (len_elem - 1)):(i*len_elem)) = str_arr(i)
          end do
        end subroutine string_array_to_string

    end subroutine broadcast

    subroutine recv(t)
        use mpi
        class(geometry) :: t
        integer :: ierr, rank
        character, dimension(:), allocatable :: atoms_string
        
        call t%dealloc_arrays()
        call mpi_comm_rank(MPI_COMM_WORLD, rank, ierr)

        call mpi_recv(t%nat, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        if (ierr /= 0) print*, 'sfadsfs'

        allocate(t%pos(3, t%nat), t%atomnames(t%nat))
        allocate(atoms_string(2*t%nat))

        call mpi_recv(t%pos, 3*t%nat, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        if (ierr /= 0) print*, 'sngjn;o'
        call mpi_recv(t%lat, 9, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr) 
        if (ierr /= 0) print*, 'xcvb i'
        call mpi_recv(atoms_string, 2*t%nat, MPI_CHARACTER, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        if (ierr /= 0) print*, 'e[g]ca;'

        call string_to_string_array(atoms_string, t%atomnames, t%nat, 2)

    contains
        subroutine string_to_string_array(str, str_arr, len_array, len_elem)
          implicit none
          integer, intent(in) :: len_array
          integer, intent(in) :: len_elem
          character, dimension(len_array*len_elem), intent(in) :: str
          character(len=len_elem), dimension(len_array), intent(out) :: str_arr
          integer :: i
          character(len=len_array*len_elem) :: tempstr
          do i = 1, len_array*len_elem, 1
            tempstr(i:i) = str(i)
          end do
          do i = 1, len_array, 1
            str_arr(i) = tempstr((i*len_elem - (len_elem - 1)):(i*len_elem))
          end do
        end subroutine string_to_string_array

    end subroutine recv
    
end module atomic_structure