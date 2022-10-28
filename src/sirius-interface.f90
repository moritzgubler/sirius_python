module mod_sirius_energy
  use sirius
  implicit none
  type(sirius_context_handler), private:: handler
  type(sirius_kpoint_set_handler), private :: kset
  type(sirius_ground_state_handler), private :: dft
  integer, private :: ntypes
  character(len=2), dimension(:), allocatable, private :: atom_types
  real(8), private :: pressure = 0.0
  logical, private :: pressure_set = .FALSE.
  logical, private :: is_setup = .FALSE.
  integer, private :: size_cluster
  integer, private :: rank
  logical, private :: set_up_first_time = .TRUE.
  integer, private, parameter :: wake_sig = 666, shutdown_sig = 777
  integer, private :: nat_sirius
  integer, private, dimension(3) :: k_points_sir
  integer, private, dimension(3) :: k_offset_sir
  logical, private :: k_sym = .true.
  character(len=2), allocatable, dimension(:) :: atomnames_sir
  character(len=200) :: json_dir_sir
  real(8), private :: pw_cutoff_sir
  real(8), private :: gk_cutoff_sir
  real(8), private :: vol_setup

  !! declare visibility of subroutines:
  public :: enthalpyforces_sirius
  public :: energyandforces_sirius
  public :: setup_sirius
  public :: set_pressure_sirius
  public :: get_pressure_sirius
  public :: set_pressure_si_sirius
  public :: exit_sirius
  public :: reset_sirius
  public :: shutdown_sirius_interface
  public :: setup_sirius_internal
  public :: set_slave_idle
  public :: recv_setup_dat_on_slave
  public :: wakeup_slave
  public :: string_array_to_string
  public :: string_to_string_array
  public :: send_setup_dat_to_slave
  public :: slave_loop
  public :: start_cal_on_slave
  public :: calc_dv_dalat
  public :: calc_vol
  public :: calc_deralat
  public :: get_atom_types
  public :: cart2frac_sirius
contains

  !! calculates enthalpy and forces. To use this subroutine this module first
  !! needs to be set up by calling the set up subroutine (with all mpi processes).
  !! All variables are in atomic units
  subroutine enthalpyforces_sirius(nat, rxyz, fxyz, enth, alat, deralat)
    implicit none
    integer :: nat
    !! Number of atoms.
    real(8), dimension(3, nat), intent(in) :: rxyz
    !! Positions of the atoms.
    real(8), dimension(3, nat) :: fxyz
    !! Forces acting on the atoms.
    real(8) :: enth
    !! enthalpy H=E+p*V of the system.
    real(8), dimension(3, 3) :: alat
    !! Lattice vectors of the periodic cell.
    real(8), dimension(3, 3) :: deralat
    !! Negative derivative of the enthalpy with respect to the lattice vectors.
    real(8) :: etot
    !! derrivative d (det(alat)) / d (alat)
    real(8), dimension(3, 3) :: dv_dalat
    real(8) :: vol
    if (.not. pressure_set) then
      print *, "pressure was not set for the enthalpy calculation"
      stop "call set_pressure in mod sirius energy"
    end if
    call energyandforces_sirius(nat, rxyz, fxyz, etot, alat, deralat)
    call calc_dv_dalat(alat, dv_dalat)
    call calc_vol(alat, vol)
    deralat = deralat - pressure*dv_dalat
    enth = etot + pressure*vol
  end subroutine enthalpyforces_sirius

  !! calculates energy and forces. To use this subroutine this module first
  !! needs to be set up by calling the set up subroutine (with all mpi processes).
  !! All variables are in atomic units
  subroutine energyandforces_sirius(nat, rxyz, fxyz, etot, alat, deralat)
    use mpi
    implicit none
    integer, intent(in) :: nat
    real(8), intent(in), dimension(3, nat) :: rxyz
    real(8), dimension(3, nat), intent(out) :: fxyz
    real(8), intent(out) :: etot
    real(8), dimension(3, 3), intent(in) :: alat
    real(8), dimension(3, 3), intent(out) :: deralat
    real(8), dimension(3, 3) :: stress
    integer :: i
    real(8), dimension(3, nat) :: xyzred
    real(8) :: scf_correction
    integer :: sir_err
    !real(8) :: vol

    !if ( rank == 0 ) then
    !  call calc_vol(alat, vol)
    !  if ( abs(vol / vol_setup - 1 ) > 0.1 ) then
    !    print*, 'automatic reset due to volume change'
    !    call reset_sirius(nat, rxyz, alat)
    !  end if
    !  if ( vol <= 0 ) then
    !    print*, "warning, negative determinant of lattice cell.", vol
    !  end if
    !end if

    if (size_cluster > 1) then !! program uses several mpi processes.
      !! invoke mpi slave processes
      if (rank == 0) then ! master process
        call start_cal_on_slave(nat, rxyz, alat)
      end if
    end if
    if (.not. is_setup) then
      stop "sirius is not setup."
    end if
    if (nat /= nat_sirius) then
      print *, "nat does not equal nat_sirius!"
      print *, "nat:", nat, "nat_sirius", nat_sirius
    end if
    call cart2frac_sirius(nat, alat, rxyz, xyzred)
    do i = 1, nat, 1
      call sirius_set_atom_position(handler, i, xyzred(:, i), error_code=sir_err)
      if(sir_err /= 0) stop "error setting atom positions"
    end do
    call sirius_set_lattice_vectors(handler, alat(:, 1), alat(:, 2), alat(:, 3), error_code=sir_err)
    if(sir_err /= 0) stop "error setting lattice vectors"

    call sirius_update_ground_state(dft, sir_err)
    if (sir_err /= 0) then
      if (rank.eq.0) then
        write(*,*)"error in sirius_update_ground_state()"
        write(*,*)"lattice vectors:"
        write(*,*)alat(:, 1)
        write(*,*)alat(:, 2)
        write(*,*)alat(:, 3)
        write(*,*)"atom coordinates:"
        do i = 1, nat, 1
          write(*,*) xyzred(:, i)
        enddo
      endif
      call flush(6)
      call mpi_barrier(MPI_COMM_WORLD, i)
      call mpi_abort(MPI_COMM_WORLD, -1, i)
    endif

    call sirius_set_parameters(handler, error_code=sir_err)
    if(sir_err /= 0) stop "error setting parameters sirius"
    call sirius_find_ground_state(dft, density_tol=1d-6, energy_tol=1d-8, error_code=sir_err)
    if(sir_err /= 0) stop "error finding ground state sirius"

    call sirius_get_forces(dft, "total", fxyz, error_code=sir_err)
    if (sir_err /= 0) stop "error getting forces"
    call sirius_get_stress_tensor(dft, "total", stress, error_code=sir_err)
    if (sir_err /= 0) stop "error getting stress tensor"
    call sirius_get_energy(dft, "total", etot, error_code=sir_err)
    if(sir_err /= 0) stop "error getting energy sirius"
    call sirius_get_energy(dft, "descf", scf_correction, error_code=sir_err)
    if(sir_err /= 0) stop "error getting scf energy correction sirius"
    etot = etot + scf_correction
    ! calculate lattice derivatives of energy from the stress tensor:
    stress = -stress
    call calc_deralat(stress, alat, deralat)
  end subroutine energyandforces_sirius

  !! use this subroutine to set up the calculations with sirius.
  !! Atomic units are used in the entire module.
  !! This subroutine needs to be called with all the mpi processes that are
  !! available. Only process 0 will return from this subroutine. The others
  !! will be kept ready until another energy calculation is required.
  !! mpi_init will be called here.
  subroutine setup_sirius(nat, rxyz, alat, atomnames_in, json_dir, pw_cutoff, gk_cutoff, k_points, k_offset, json_string)
    use mpi
    implicit none
    integer, intent(in) :: nat
    real(8), intent(in) :: rxyz(3, nat)
    real(8), dimension(3, 3), intent(in) :: alat
    character, dimension(2,nat), intent(in) :: atomnames_in
    character(len=100), intent(in) :: json_dir
    real(8), intent(in) :: pw_cutoff
    !! cutoff in rydberg
    real(8), intent(in) :: gk_cutoff
    ! cutoff in rydberg
    integer, dimension(3), intent(in) :: k_points
    integer, dimension(3), intent(in) :: k_offset
    character(len=*) :: json_string
    integer :: mpi_err
    integer :: t1, t2, i

    character(len=2), dimension(nat) :: atomnames
    character :: temp_char(2)


    do i = 1, nat
      temp_char = atomnames_in(:, i)
      atomnames(i)(1:1) = temp_char(1)
      atomnames(i)(2:2) = temp_char(2)
    end do

    if ( is_setup ) then ! should not be set up already. stopping.
      print*, "setup_sirius called while being setup"
      call mpi_abort(MPI_COMM_WORLD, 1, mpi_err)
    end if
    if (set_up_first_time) then
      call mpi_init(mpi_err)
      call mpi_comm_size(MPI_COMM_WORLD, size_cluster, mpi_err)
      if (mpi_err /= 0) then
        stop "mpi_comm_size fails in set_up_sirius"
      end if
      call mpi_comm_rank(MPI_COMM_WORLD, rank, mpi_err)
      if (mpi_err /= 0) then
        stop "error getting mpi rank in setup_sirius."
      end if
      ! initialize the library; this should be called only once during the mpi execution
      call sirius_initialize(call_mpi_init=.false., error_code=mpi_err)
      if ( mpi_err /= 0 ) then
        stop "error initializing sirius"
      end if
      !set_up_first_time = .FALSE.
    end if
    if (size_cluster > 1) then ! mpi parallelizazion present
      if (rank == 0) then ! is on master process
        if (.not. set_up_first_time) then !! mpi threads are idle in set slave idle
          call wakeup_slave(nat, rxyz, alat, atomnames, json_dir, pw_cutoff, gk_cutoff, k_points, k_offset)
        end if
      end if
    end if
    set_up_first_time = .FALSE.
    call setup_sirius_internal(nat, rxyz, alat, atomnames, json_dir, pw_cutoff, gk_cutoff, k_points, k_offset)
    if (rank > 0) then ! process is slave
      !open(6,file="standard_out_garbage.out")
      call slave_loop
      stop "slave processes shouldnt be allowed leaving this module."
    end if
  end subroutine setup_sirius

  subroutine setup_sirius_internal(nat, rxyz, alat, atomnames, json_dir, pw_cutoff, gk_cutoff, k_points, k_offset, json_string)
    use mpi
    implicit none
    integer, intent(in) :: nat
    real(8), intent(in) :: rxyz(3, nat)
    real(8), intent(in) :: alat(3, 3)
    character(len=2), intent(in) :: atomnames(nat)
    character(len=*), intent(in) :: json_dir
    real(8),intent(in) :: pw_cutoff
    real(8), intent(in) :: gk_cutoff
    integer, intent(in) :: k_points(3)
    integer, intent(in) :: k_offset(3)
    character(len=*), intent(in) :: json_string
    integer :: i, nkpt
    real(8), dimension(3, nat) :: xyzred
    integer :: sir_err
    character(len=150) :: pot_file
    logical :: pot_file_exists, first_try = .TRUE.

    if (is_setup) then !! when sirius is already setup, cancel
      stop "is already set up"
    end if
    nat_sirius = nat
    k_points_sir = k_points
    k_offset_sir = k_offset
    if ( allocated(atomnames_sir) ) deallocate(atomnames_sir)
    allocate(atomnames_sir(nat))
    atomnames_sir = atomnames
    json_dir_sir = json_dir
    pw_cutoff_sir = sqrt(pw_cutoff)
    gk_cutoff_sir = sqrt(gk_cutoff)
    1234 continue
    ! create simulation context using a specified communicator
    call sirius_create_context(MPI_COMM_WORLD, handler, error_code=sir_err)
    if (sir_err /= 0) stop "error in creating context"
    !call sirius_import_parameters(handler, &
    !  '{"parameters" : {"electronic_structure_method" : "pseudopotential"},&
    !  "control" : {"verbosity" : 0, "verification" : 0}}', error_code=sir_err)
    ! parameters from anton
    call sirius_import_parameters(handler, &
    '{"parameters" : {"electronic_structure_method" : "pseudopotential"},&
      &    "control" : {"verbosity" : 0, "verification" : 0},&
      &    "mixer" : {"beta" : 0.5, "max_history" : 8, "use_hartree" : true}}', error_code=sir_err)
    if (sir_err /= 0) stop "error in sirius_import_parameters"
    ! atomic units are used everywhere
    ! plane-wave cutoffs are provided in a.u.^-1
    call sirius_set_parameters(handler, pw_cutoff=pw_cutoff_sir, gk_cutoff=gk_cutoff_sir &
                               , use_symmetry=k_sym, error_code=sir_err)
    if (sir_err /= 0) stop "error in set parameters"
    call get_atom_types(nat, atomnames)
    do i = 1, ntypes, 1
      pot_file = trim(json_dir)//trim(atom_types(i))//".json"
      inquire (file=pot_file, exist=pot_file_exists)
      if (.not. pot_file_exists) then
        print *, 'pot_file "'//trim(pot_file)//'" does not exist.'
        print*, 'rank', rank
        stop "pot file does not exist."
      end if
      call sirius_add_atom_type(handler, atom_types(i) &
                                , fname=trim(pot_file), error_code=sir_err)
      if (sir_err /= 0) then
        print *, "error in add atom_type ", i, atom_types(i)
        stop "error adding atom type"
      end if
    end do
    call sirius_set_lattice_vectors(handler, alat(:, 1), alat(:, 2), alat(:, 3))
    call cart2frac_sirius(nat, alat, rxyz, xyzred)
    do i = 1, nat, 1
      call sirius_add_atom(handler, atomnames(i), xyzred(:, i))
    end do
    ! exchange and correlation parts
    !call sirius_add_xc_functional(handler, "XC_LDA_X")
    !call sirius_add_xc_functional(handler, "XC_LDA_C_VWN")

    ! PBE functinals
    call sirius_add_xc_functional(handler, "XC_GGA_X_PBE")
    call sirius_add_xc_functional(handler, 'XC_GGA_C_PBE')

    ! initialize the simulation handler
    call sirius_initialize_context(handler, error_code=sir_err)
    if (sir_err /= 0) then
      !call save_me(nat, alat, rxyz)
      print*, "failing rank", rank, first_try
      call sleep(1)
      if ( .not. first_try  ) STOP "error calling sirius_initialize_context()"
      first_try = .FALSE.
      call sirius_free_handler(handler)
      call sleep(1)
      print*, "second time", sir_err
      goto 1234
      STOP "error calling sirius_initialize_context()"
    end if
    call sirius_create_kset_from_grid(handler, k_points_sir, k_offset_sir, k_sym, kset, error_code=sir_err)
    if (sir_err /= 0) stop "error creating k_grid"

    call sirius_create_ground_state(kset, dft, error_code=sir_err)
    if (sir_err /= 0) stop "error creating ground state"
    call sirius_get_num_kpoints(kset, nkpt, error_code=sir_err)
    if (sir_err /= 0) stop "error getting number of k-points"
    if (rank.eq.0) then
      write(*,*)'number of k-points : ', nkpt
    endif

    !! set volume of cell during initialization
    call calc_vol(alat, vol_setup)
    is_setup = .TRUE.
  end subroutine setup_sirius_internal

  !! slave precesses will wait in this subroutine after library is shut down
  subroutine set_slave_idle
    use mpi
    implicit none
    integer :: wakeup
    integer :: mpi_err
    if (rank == 0) then
      print *, "set slave idle was called with master process"
      call mpi_abort(MPI_COMM_WORLD, 1, mpi_err)
    end if
    call mpi_recv(wakeup, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpi_err)
    if (wakeup == wake_sig) then !! get setup data and return
      call recv_setup_dat_on_slave
      return
    end if
    if (wakeup == shutdown_sig) then
      call mpi_finalize(mpi_err)
      stop
    end if
    print *, "received wrong wakup signal: ", wakeup
    call mpi_abort(MPI_COMM_WORLD, 1, mpi_err)
  end subroutine set_slave_idle

  subroutine recv_setup_dat_on_slave
    use mpi
    implicit none
    integer :: nat
    real(8), dimension(:), allocatable :: rxyz
    real(8), dimension(3, 3) :: alat
    character, dimension(:), allocatable :: atomnames
    character(len=2), dimension(:), allocatable :: atomstring
    integer :: len_json
    character*200 :: json_dir
    real(8) :: pw_cutoff
    real(8) :: gk_cutoff
    integer, dimension(3) :: k_points
    integer, dimension(3) :: k_offset
    integer :: mpi_err

    if (rank == 0) then
      print *, "recv setup_dat_on_slave was called from master process."
      call mpi_abort(MPI_COMM_WORLD, 1, mpi_err)
    end if
    if (is_setup) then
      print *, "recv_setup_dat_on_slave was called while mod sirius was set up."
      call mpi_abort(MPI_COMM_WORLD, 1, mpi_err)
    end if
    if (set_up_first_time) then
      print *, "recv_setup_dat_on_slave was called before it was setup first time"
      call mpi_abort(MPI_COMM_WORLD, 1, mpi_err)
    end if

    call mpi_recv(nat, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpi_err)
    if (mpi_err /= 0) then
      print *, "error receiving nat in recv_setup_dat_on_slave."
      call mpi_abort(MPI_COMM_WORLD, 1, mpi_err)
    end if
    allocate (rxyz(3*nat))
    allocate (atomnames(2*nat))
    allocate (atomstring(nat))
    call mpi_recv(rxyz, 3*nat, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpi_err)
    if (mpi_err /= 0) then
      print *, "error receiving rxyz in recv_setup_dat_on_slave."
      call mpi_abort(MPI_COMM_WORLD, 1, mpi_err)
    end if
    call mpi_recv(alat, 9, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpi_err)
    if (mpi_err /= 0) then
      print *, "error receiving alat in recv_setup_dat_on_slave."
      call mpi_abort(MPI_COMM_WORLD, 1, mpi_err)
    end if
    call mpi_recv(atomnames, 2*nat, MPI_CHARACTER, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpi_err)
    if (mpi_err /= 0) then
      print *, "error receiving atomnames in recv_setup_dat_on_slave."
      call mpi_abort(MPI_COMM_WORLD, 1, mpi_err)
    end if
    call mpi_recv(len_json, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpi_err)
    if (mpi_err /= 0) then
      print *, "error receiving len_json in recv_setup_dat_on_slave."
      call mpi_abort(MPI_COMM_WORLD, 1, mpi_err)
    end if
    if (len_json > 200) then
      print *, "length of jsondir string cant be longer than 200 characters. aborting"
      call mpi_abort(MPI_COMM_WORLD, 1, mpi_err)
    end if
    json_dir = ""
    call mpi_recv(json_dir(1:len_json), len_json, MPI_CHARACTER, 0, 1 &
                  , MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpi_err)
    if (mpi_err /= 0) then
      print *, "error receiving jsondir in recv_setup_dat_on_slave."
      call mpi_abort(MPI_COMM_WORLD, 1, mpi_err)
    end if
    call mpi_recv(pw_cutoff, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpi_err)
    if (mpi_err /= 0) then
      print *, "error receiving pw_cutoff in recv_setup_dat_on_slave."
      call mpi_abort(MPI_COMM_WORLD, 1, mpi_err)
    end if
    call mpi_recv(gk_cutoff, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpi_err)
    if (mpi_err /= 0) then
      print *, "error receiving gk_cutoff in recv_setup_dat_on_slave."
      call mpi_abort(MPI_COMM_WORLD, 1, mpi_err)
    end if
    call mpi_recv(k_points, 3, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpi_err)
    if (mpi_err /= 0) then
      print *, "error receiving k_points in recv_setup_dat_on_slave."
      call mpi_abort(MPI_COMM_WORLD, 1, mpi_err)
    end if
    call mpi_recv(k_offset, 3, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpi_err)
    if (mpi_err /= 0) then
      print *, "error receiving k_offset in recv_setup_dat_on_slave."
      call mpi_abort(MPI_COMM_WORLD, 1, mpi_err)
    end if
    call string_to_string_array(atomnames, atomstring, nat, 2)
    call setup_sirius_internal(nat, rxyz, alat, atomnames, json_dir, pw_cutoff, gk_cutoff, k_points, k_offset)
  end subroutine recv_setup_dat_on_slave

  !! master process can use this subroutine to get idle slave threads to call
  !! setup_sirius.
  subroutine wakeup_slave(nat, rxyz, alat, atomnames, json_dir, pw_cutoff, gk_cutoff, k_points, k_offset)
    use mpi
    implicit none
    integer, intent(in) :: nat
    real(8), dimension(3, nat), intent(in) :: rxyz
    real(8), dimension(3, 3), intent(in) :: alat
    character(len=2), dimension(nat), intent(in) :: atomnames
    character(len=*) :: json_dir
    real(8) :: pw_cutoff
    real(8) :: gk_cutoff
    integer, dimension(3), intent(in) :: k_points
    integer, dimension(3), intent(in) :: k_offset
    integer :: i
    integer :: mpi_err
    character(len=2*nat) :: atomstring
    if (rank /= 0) then
      print *, "wake up slave from slave thread called."
      call mpi_abort(MPI_COMM_WORLD, 1, mpi_err)
    end if
    do i = 1, size_cluster - 1, 1
      call mpi_send(wake_sig, 1, MPI_INT, i, 1, MPI_COMM_WORLD, mpi_err)
      if (mpi_err /= 0) then
        print *, "error sending wakeup signal to slave with rank ", i
        call mpi_abort(MPI_COMM_WORLD, 1, mpi_err)
      end if
    end do
    call string_array_to_string(atomstring, atomnames, nat, 2)
    call send_setup_dat_to_slave(nat, rxyz, alat, atomstring, trim(json_dir) &
                                 , len_trim(json_dir), pw_cutoff, gk_cutoff, k_points, k_offset)
  end subroutine wakeup_slave

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

  subroutine send_setup_dat_to_slave(nat, rxyz, alat, atomnames, json_dir, len_json &
                                     , pw_cutoff, gk_cutoff, k_points, k_offset)
    use mpi
    implicit none
    integer, intent(in) :: nat
    real(8), dimension(3, nat), intent(in) :: rxyz
    real(8), dimension(3, 3), intent(in) :: alat
    character(len=2*nat), intent(in) :: atomnames
    integer :: len_json
    character(len=len_json) :: json_dir
    real(8) :: pw_cutoff
    real(8) :: gk_cutoff
    integer, dimension(3), intent(in) :: k_points
    integer, dimension(3), intent(in) :: k_offset
    integer :: i
    integer :: mpi_err
    do i = 1, size_cluster - 1, 1
      call mpi_send(nat, 1, MPI_INT, i, 1, MPI_COMM_WORLD, mpi_err)
      call mpi_send(rxyz, 3*nat, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, mpi_err)
      call mpi_send(alat, 9, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, mpi_err)
      call mpi_send(atomnames, 2*nat, MPI_CHARACTER, i, 1, MPI_COMM_WORLD, mpi_err)
      call mpi_send(len_json, 1, MPI_INT, i, 1, MPI_COMM_WORLD, mpi_err)
      call mpi_send(json_dir, len_json, MPI_CHARACTER, i, 1, MPI_COMM_WORLD, mpi_err)
      call mpi_send(pw_cutoff, 1, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, mpi_err)
      call mpi_send(gk_cutoff, 1, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, mpi_err)
      call mpi_send(k_points, 3, MPI_INT, i, 1, MPI_COMM_WORLD, mpi_err)
      call mpi_send(k_offset, 3, MPI_INT, i, 1, MPI_COMM_WORLD, mpi_err)
    end do
  end subroutine send_setup_dat_to_slave

  !! this subroutine will send a shutdown signal to all slave processes which will terminate them.
  !! master thread will call mpi_finalize. the master process is still available
  !! after this subroutine.
  subroutine shutdown_sirius_interface
    use mpi
    implicit none
    integer :: i
    integer :: mpi_err
    if (is_setup) then
      print *, "shutdown sirius was called befor sirius was exited."
      print *, "this forbidden action is cancelled and program will continue in a safe way."
      return
    end if
    if (rank /= 0) then
      print *, "shutdown sirius interface called from slave thread. this is forbidden."
      call mpi_abort(MPI_COMM_WORLD, 1, i)
    end if
    if (set_up_first_time) then
      print *, "shutdown sirius interface was called before sirius interface was initialized first time"
      print *, " this call is being ignored and the program will continue"
      return
    end if
    do i = 1, size_cluster - 1, 1
      call mpi_send(shutdown_sig, 1, MPI_INT, i, 1, MPI_COMM_WORLD, mpi_err)
    end do
    ! from the MPI documentaion:
    ! All processes must call this routine before exiting.
    call mpi_finalize(i)
    if (i /= 0) then
      print *, "mpi_finalize failed on master with error code: ", i
    end if
  end subroutine shutdown_sirius_interface

  subroutine slave_loop!(nat)
    use mpi
    implicit none
    !integer :: nat
    real(8) :: etot
    real(8), dimension(3, nat_sirius) :: rxyz
    real(8), dimension(3, nat_sirius) :: fxyz
    real(8), dimension(3, 3) :: alat
    real(8), dimension(3, 3) :: deralat
    integer :: do_ener_cal
    integer :: mpi_err
    if (rank == 0) then
      stop "slave loop was called on master process"
    end if
    do
      !! if shut down, stop thread.
      call mpi_recv(do_ener_cal, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpi_err)
      if (mpi_err /= 0) then
        stop "error receiving enercal tag in slave"
      end if
      if (do_ener_cal == 0) then ! no more energy calculations needed. shutting down
        call exit_sirius()
        call set_slave_idle()
        cycle
      end if
      if (do_ener_cal /= 1) then
        print *, "enercal tag", do_ener_cal
        stop "do_enercal tag not equal to 0 or 1 in slave"
      end if
      !! receive data from master thread.
      !! receive rxyz and alat
      call mpi_recv(rxyz, 3*nat_sirius, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpi_err)
      if (mpi_err /= 0) then
        stop "error receiving positions in slave"
      end if
      call mpi_recv(alat, 9, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpi_err)
      if (mpi_err /= 0) then
        stop "error receiving lattice vectors in slave"
      end if
      !! do energy calculation.
      call energyandforces_sirius(nat_sirius, rxyz, fxyz, etot, alat, deralat)
    end do
  end subroutine slave_loop

  subroutine start_cal_on_slave(nat, rxyz, alat)
    use mpi
    implicit none
    integer :: nat
    real(8), dimension(3, nat) :: rxyz
    real(8), dimension(3, 3) :: alat
    integer :: mpi_err
    integer :: i
    if (size_cluster == 1) then
      return
    end if
    if (rank /= 0) then
      stop "start_cal_on_slave was executed by a slave process"
    end if
    do i = 1, size_cluster - 1, 1
      !! send 1 (true) to slave. this signals the slave to get ready for an energy calculation.
      !! sening a 0 (false) would trigger the slave to shut down sirius library.
      call mpi_send(1, 1, MPI_INT, i, 1, MPI_COMM_WORLD, mpi_err)
      if (mpi_err /= 0) then
        stop "error sending force eval tag to slave."
      end if
      call MPI_SEND(rxyz, 3*nat, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, mpi_err)
      if (mpi_err /= 0) then
        stop "error sending positions to slave"
      end if
      call MPI_SEND(alat, 9, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, mpi_err)
      if (mpi_err /= 0) then
        stop "error sending lattice vectors to slave"
      end if
    end do
  end subroutine start_cal_on_slave

  !! Calculates the derivative of the volume of the cell with respect
  !! to the lattice vectors.
  subroutine calc_dv_dalat(alat, dv_dalat)
    implicit none
    !! Lattice vectors of the periodic cell.
    real(8), dimension(3, 3), intent(in) :: alat
    !! Derivative of the volume of the cell (determinant) with respect to the
    !! lattice parameters.
    real(8), dimension(3, 3), intent(out) :: dv_dalat
    dv_dalat(1, 1) = alat(2, 2)*alat(3, 3) - alat(2, 3)*alat(3, 2)
    dv_dalat(2, 1) = alat(1, 3)*alat(3, 2) - alat(1, 2)*alat(3, 3)
    dv_dalat(3, 1) = alat(1, 2)*alat(2, 3) - alat(2, 2)*alat(1, 3)
    dv_dalat(1, 2) = alat(2, 3)*alat(3, 1) - alat(2, 1)*alat(3, 3)
    dv_dalat(2, 2) = alat(1, 1)*alat(3, 3) - alat(3, 1)*alat(1, 3)
    dv_dalat(3, 2) = alat(1, 3)*alat(2, 1) - alat(1, 1)*alat(2, 3)
    dv_dalat(1, 3) = alat(2, 1)*alat(3, 2) - alat(3, 1)*alat(2, 2)
    dv_dalat(2, 3) = alat(1, 2)*alat(3, 1) - alat(1, 1)*alat(3, 2)
    dv_dalat(3, 3) = alat(1, 1)*alat(2, 2) - alat(2, 1)*alat(1, 2)
  end subroutine calc_dv_dalat

  subroutine calc_vol(alat, vol)
    implicit none
    real(8), dimension(3, 3), intent(in) :: alat
    real(8), intent(out) :: vol
    vol = alat(1, 1)*alat(2, 2)*alat(3, 3) - alat(1, 1)*alat(2, 3)*alat(3, 2) - &
          alat(1, 2)*alat(2, 1)*alat(3, 3) + alat(1, 2)*alat(2, 3)*alat(3, 1) + &
          alat(1, 3)*alat(2, 1)*alat(3, 2) - alat(1, 3)*alat(2, 2)*alat(3, 1)
  end subroutine calc_vol

  subroutine calc_deralat(stress, alat, deralat)
    implicit none
    real(8), dimension(3, 3), intent(in) :: stress
    real(8), dimension(3, 3), intent(in) :: alat
    real(8), dimension(3, 3), intent(out) :: deralat
    real(8), dimension(3, 3) :: alatinv
    real(8) :: vol
    real(8) :: div
    call calc_vol(alat, vol)
    div = 1.d0/vol
    alatinv(1, 1) = (alat(2, 2)*alat(3, 3) - alat(2, 3)*alat(3, 2))*div
    alatinv(1, 2) = -(alat(1, 2)*alat(3, 3) - alat(1, 3)*alat(3, 2))*div
    alatinv(1, 3) = (alat(1, 2)*alat(2, 3) - alat(1, 3)*alat(2, 2))*div
    alatinv(2, 1) = -(alat(2, 1)*alat(3, 3) - alat(2, 3)*alat(3, 1))*div
    alatinv(2, 2) = (alat(1, 1)*alat(3, 3) - alat(1, 3)*alat(3, 1))*div
    alatinv(2, 3) = -(alat(1, 1)*alat(2, 3) - alat(1, 3)*alat(2, 1))*div
    alatinv(3, 1) = (alat(2, 1)*alat(3, 2) - alat(2, 2)*alat(3, 1))*div
    alatinv(3, 2) = -(alat(1, 1)*alat(3, 2) - alat(1, 2)*alat(3, 1))*div
    alatinv(3, 3) = (alat(1, 1)*alat(2, 2) - alat(1, 2)*alat(2, 1))*div
    alatinv = transpose(alatinv)
    deralat = vol*matmul(stress, alatinv)
  end subroutine calc_deralat

  !! delete all variables that have been setup. after execution sirius will need
  !! to be set up again. The mpi processes are being held in this module.
  subroutine exit_sirius()
    use mpi
    implicit none
    !logical :: call_mpi_fin = .FALSE.
    !! if FALSE, mpi_finalize must be called after shutdown
    integer :: i
    integer :: mpi_err
    if ( .not. is_setup ) then
      !print*, "exit sirius called without sirius being initialized. leaving subroutine again."
      return
    end if
    if (size_cluster > 1) then
      if (rank == 0) then ! master process
        do i = 1, size_cluster - 1, 1
          call mpi_send(0, 1, MPI_INT, i, 1, MPI_COMM_WORLD, mpi_err)
          if (mpi_err /= 0) then
            stop "error sending abort signal to slave thread"
          end if
        end do
      end if
    end if
    nat_sirius = -1
    call sirius_free_handler(dft)
    call sirius_free_handler(kset)
    call sirius_free_handler(handler)
    !call sirius_finalize(call_mpi_fin=call_mpi_fin)
    !if (rank.eq.0) then
    !  call sirius_print_timers(flatten=.false.)
    !endif
    pressure = 0.0
    is_setup = .FALSE.
    pressure_set = .FALSE.
  end subroutine exit_sirius

  !! delete all variables that have been setup and sets it up again with the previous parameters
  !! but updated postions and lattice vectors.
  subroutine reset_sirius(nat, rxyz, alat)
    use mpi
    implicit none
    integer, intent(in) :: nat
    real(8), intent(in) :: rxyz(3,nat), alat(3, 3)
    character(len=2), dimension(nat) :: symb
    real(8)  :: p_priv
    real(8) :: pw_cutoff_rydberg, gk_cutoff_rydberg
    if ( .not. is_setup ) then
      !print*, "exit sirius called without sirius being initialized. leaving subroutine again."
      return
    end if
    if (size_cluster > 1) then
      if (rank /= 0) then ! master process
        stop "reset_sirius must be called with master process"
      end if
    end if
    call exit_sirius
    print*, "resetting"!, nat, atomnames_sir
    symb = atomnames_sir
    pw_cutoff_rydberg = pw_cutoff_sir**2
    gk_cutoff_rydberg = gk_cutoff_sir**2
    print*, 'cutoffs, pw, gk', pw_cutoff_rydberg, gk_cutoff_rydberg
    call setup_sirius(nat, rxyz, alat, symb, json_dir_sir, pw_cutoff_rydberg&
      , gk_cutoff_rydberg, k_points_sir, k_offset_sir)
    p_priv = pressure
    call set_pressure_sirius(p_priv)
  end subroutine reset_sirius

  subroutine get_atom_types(nat, atomnames)
    integer, intent(in) :: nat
    character(len=2), dimension(nat), intent(in) :: atomnames
    integer :: i, j
    character(len=2), dimension(nat) :: an_copy

    ntypes = 0
    an_copy = atomnames
    do i = 1, nat, 1
      if (len_trim(an_copy(i)) /= 0) then ! new atom species found. delete all the following from list
        ntypes = ntypes + 1
        do j = i + 1, nat, 1
          if (an_copy(j) == an_copy(i)) then
            an_copy(j) = ""
          end if
        end do
      end if
    end do
    if (allocated(atom_types)) then
      deallocate (atom_types)
    end if
    allocate (atom_types(ntypes))
    j = 0
    do i = 1, nat, 1
      if (an_copy(i) /= "") then
        j = j + 1
        atom_types(j) = an_copy(i)
      end if
    end do
  end subroutine get_atom_types

  !! converts cartesian coordinates rxyz to reduced coordinates xyzred
  subroutine cart2frac_sirius(nat, alat, rxyz, xyzred)
    implicit none
    !! Number of Atoms
    integer, intent(in) :: nat
    !! Lattice Vectors.
    real(8), intent(in), dimension(3, 3) :: alat
    !! Position of the Atoms in cartesian coorinates.
    real(8), intent(in), dimension(3, nat) :: rxyz
    !! Position of the Atoms in reduced coordinates.
    real(8), intent(out), dimension(3, nat) :: xyzred

    !private variables
    integer :: iat, l
    real(8), dimension(3, 3) :: alatinv
    real(8) :: div

    div = alat(1, 1)*alat(2, 2)*alat(3, 3) - alat(1, 1)*alat(2, 3)*alat(3, 2) - &
          alat(1, 2)*alat(2, 1)*alat(3, 3) + alat(1, 2)*alat(2, 3)*alat(3, 1) + &
          alat(1, 3)*alat(2, 1)*alat(3, 2) - alat(1, 3)*alat(2, 2)*alat(3, 1)
    div = 1.d0/div
    alatinv(1, 1) = (alat(2, 2)*alat(3, 3) - alat(2, 3)*alat(3, 2))*div
    alatinv(1, 2) = -(alat(1, 2)*alat(3, 3) - alat(1, 3)*alat(3, 2))*div
    alatinv(1, 3) = (alat(1, 2)*alat(2, 3) - alat(1, 3)*alat(2, 2))*div
    alatinv(2, 1) = -(alat(2, 1)*alat(3, 3) - alat(2, 3)*alat(3, 1))*div
    alatinv(2, 2) = (alat(1, 1)*alat(3, 3) - alat(1, 3)*alat(3, 1))*div
    alatinv(2, 3) = -(alat(1, 1)*alat(2, 3) - alat(1, 3)*alat(2, 1))*div
    alatinv(3, 1) = (alat(2, 1)*alat(3, 2) - alat(2, 2)*alat(3, 1))*div
    alatinv(3, 2) = -(alat(1, 1)*alat(3, 2) - alat(1, 2)*alat(3, 1))*div
    alatinv(3, 3) = (alat(1, 1)*alat(2, 2) - alat(1, 2)*alat(2, 1))*div

    do iat = 1, nat
      xyzred(1, iat) = alatinv(1, 1)*rxyz(1, iat) + alatinv(1, 2)*rxyz(2, iat) + alatinv(1, 3)*rxyz(3, iat)
      xyzred(2, iat) = alatinv(2, 1)*rxyz(1, iat) + alatinv(2, 2)*rxyz(2, iat) + alatinv(2, 3)*rxyz(3, iat)
      xyzred(3, iat) = alatinv(3, 1)*rxyz(1, iat) + alatinv(3, 2)*rxyz(2, iat) + alatinv(3, 3)*rxyz(3, iat)
    end do
    do iat = 1, nat
      do l = 1, 3
        xyzred(l, iat) = modulo(xyzred(l, iat), 1.d0)
      end do
    end do
  end subroutine cart2frac_sirius

  !! set pressure for enthalpy calc. p_set is in hartree.
  subroutine set_pressure_sirius(p_set)
    implicit none
    real(8), intent(in) :: p_set
    pressure_set = .TRUE.
    pressure = p_set
  end subroutine set_pressure_sirius

  real(8) function get_pressure_sirius()
    get_pressure_sirius = pressure
  end function get_pressure_sirius

  subroutine set_pressure_si_sirius(p_si)
    implicit none
    real(8), intent(in) :: p_si
    real(8) :: p_convert = (1.0/2.9421015697)*10.0**(-13)
    call set_pressure_sirius(p_si*p_convert)
  end subroutine set_pressure_si_sirius

  subroutine rand_rot(nat, rxyz, alat)
    !! performs a random rotation of rxyz and alat.
    integer, intent(in) :: nat
    real(8) :: rxyz(3,nat)
    real(8) :: alat(3,3)
    real(8) :: q(3,3)
    integer, parameter :: n = 3
    integer :: i, j
    integer, parameter :: lwork = 1024
    real(8) :: work(lwork), tau(n)
    integer :: info

    do i = 1, 3, 1
      do j = 1, 3, 1
        call random_number(q(i,j))
      end do
    end do
    call dgeqrf(n, n, q, n, tau, work, lwork, info)
    call dorgqr(n, n, n, q, n, tau, work, lwork, info)
    alat = matmul(q, alat)
    do i = 1, nat, 1
      rxyz(:,i) = matmul(q, rxyz(:,i))
    end do
  end subroutine rand_rot

end module mod_sirius_energy
