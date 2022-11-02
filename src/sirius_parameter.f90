module sirius_parameters
    use atomic_structure
    implicit none
    
    type sirius_params
        type(geometry) :: struc
        integer :: ntype
        character(len=2), dimension(:), allocatable :: at_types
        real(8) :: pw_cutoff
        real(8) :: gk_cutoff
        integer :: k_grid(3)
        integer :: k_shift(3)
        character(len=250), dimension(:), allocatable :: pseudo_potentials
        character(len=100), dimension(2) :: xc_functionals
        character(len=8000) :: json_string

        contains

        procedure :: initialize_sirius_params
        procedure :: deallocate_sirius_params
        procedure :: broadcast_sirius_params
        procedure :: recv_sirius_params

    end type sirius_params

contains
    subroutine initialize_sirius_params(t, nat, ntype, at_types, pos, lat, atomnames, pw_cutoff&
        , gk_cutoff, k_grid, k_shift, pseudo_potentials, xc_functionals, json_string)
        class(sirius_params) :: t
        integer :: nat
        real(8), dimension(3, nat) :: pos
        real(8), dimension(3, 3) :: lat
        character(len=2), dimension(nat) :: atomnames
        integer :: ntype
        character(len=2), dimension(ntype) :: at_types
        real(8) :: pw_cutoff
        real(8) :: gk_cutoff
        integer :: k_grid(3)
        integer :: k_shift(3)
        character(len=250), dimension(ntype) :: pseudo_potentials
        character(len=100), dimension(2) :: xc_functionals
        character(len=8000) :: json_string

        call t%struc%initialize(nat, pos, lat, atomnames)

        if ( allocated(t%at_types) ) deallocate(t%at_types)
        if ( allocated(t%pseudo_potentials)) deallocate(t%pseudo_potentials)
        allocate(t%at_types(ntype))
        allocate(t%pseudo_potentials(ntype))

        t%ntype = ntype
        t%at_types = at_types
        t%pw_cutoff = pw_cutoff
        t%gk_cutoff = gk_cutoff
        t%k_grid = k_grid
        t%k_shift = k_shift
        t%pseudo_potentials = pseudo_potentials
        t%xc_functionals = xc_functionals
        t%json_string = json_string

    end subroutine initialize_sirius_params

    subroutine deallocate_sirius_params(t)
        class(sirius_params) :: t

        if ( allocated(t%at_types) ) deallocate(t%at_types)
        if ( allocated(t%pseudo_potentials) ) deallocate(t%pseudo_potentials)

    end subroutine deallocate_sirius_params

    subroutine broadcast_sirius_params(t)
        use mpi
        class(sirius_params) :: t
        integer :: ierr, rank
        character(len=2*t%ntype) :: types_string
        character(len=t%ntype * 250) :: pseudo_string
        character(len=2*100) :: xc_string

        call mpi_comm_rank(MPI_COMM_WORLD, rank, ierr)
        
        call string_array_to_string(types_string, t%at_types, t%ntype, 2)
        call string_array_to_string(pseudo_string, t%pseudo_potentials, t%ntype, 250)
        call string_array_to_string(xc_string, t%xc_functionals, 2, 100)

        call mpi_bcast(t%ntype, 1, MPI_INT, rank, MPI_COMM_WORLD,ierr)
        if (ierr /= 0) stop './,kolp;'
        call mpi_bcast(types_string, t%ntype * 2, MPI_CHARACTER, rank, MPI_COMM_WORLD, ierr)
        if (ierr /= 0 ) stop 'kl;q4tjkln'
        call mpi_bcast(t%pw_cutoff, 1, MPI_DOUBLE, rank,MPI_COMM_WORLD, ierr)
        if (ierr /= 0) stop '/.opkma][p'
        call mpi_bcast(t%gk_cutoff, 1, MPI_DOUBLE, rank, MPI_COMM_WORLD, ierr)
        if (ierr /= 0 ) stop ';lo49wem'
        call mpi_bcast(t%k_grid, 3, MPI_INT, rank, MPI_COMM_WORLD, ierr)
        if (ierr /= 0) stop '45t6yhn weoru8ifaj'
        call mpi_bcast(t%k_shift, 3, MPI_INT, rank, MPI_COMM_WORLD, ierr)
        if (ierr /= 0) stop '45ta;sidjf6yhn weoru8ifaj'
        call mpi_bcast(pseudo_string, t%ntype*250, MPI_CHARACTER, rank, MPI_COMM_WORLD, ierr)
        if (ierr /= 0 ) stop ';p96uikjnqa'
        call mpi_bcast(xc_string, 2*100, MPI_CHARACTER, rank, MPI_COMM_WORLD, ierr)
        if(ierr /= 0) stop '005klc'
        call mpi_bcast(t%json_string, 8000, MPI_CHARACTER, rank, MPI_CHARACTER, ierr)
        if (ierr /= 0 ) stop 'poiger4235t'
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
    end subroutine broadcast_sirius_params

    subroutine recv_sirius_params(t)
        use mpi
        class(sirius_params) :: t
        integer :: ierr, rank
        character, allocatable, dimension(:) :: types_string, pseudo_string
        character(len=2*100) :: xc_string

        call t%deallocate_sirius_params()

        call mpi_comm_rank(MPI_COMM_WORLD, rank, ierr)

        call mpi_recv(t%ntype, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        if (ierr /= 0) stop 'jklajkl'

        allocate(types_string(2*t%ntype))
        allocate(pseudo_string(t%ntype * 250))

        call mpi_recv(types_string, 2*t%ntype, MPI_CHARACTER, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        if (ierr /= 0) stop 'agrfjldks;alsfgkdj;'
        call string_to_string_array(types_string, t%at_types, t%ntype, 2)

        call mpi_recv(t%pw_cutoff, 1, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        if (ierr/= 0) stop '[]p0-jkn'
        call mpi_recv(t%gk_cutoff, 1, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        if (ierr/= 0) stop '[]yarafp0-jkn'
        call mpi_recv(t%k_grid, 3, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        if (ierr/=0) stop ';m,asfd1'
        call mpi_recv(t%k_shift, 3, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        if (ierr/=0) stop ';m,asfoodf7458-2'

        call mpi_recv(pseudo_string, t%ntype*250, MPI_CHARACTER, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        if (ierr /= 0 ) stop ';o523497867'
        call string_to_string_array(pseudo_string, t%pseudo_potentials, t%ntype, 250)

        call mpi_recv(xc_string, 2*100, MPI_CHARACTER, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        if (ierr /= 0 ) stop '0wersdfva'
        call string_to_string_array(xc_string, t%xc_functionals, 100, 2)

        call mpi_recv(t%json_string, 8000, MPI_CHARACTER, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        if (ierr /= 0) stop 'error getting json'

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
    end subroutine recv_sirius_params

end module sirius_parameters