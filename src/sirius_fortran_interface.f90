module sirius_interface
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




contains

    subroutine setup_sirius(nat, ntype, types, pos, lat, atom_names, pw_cutoff, gk_cutoff, k_grid,&
         k_shift, pseudo_potentials, xc_functionals, json_string)
        
        
    end subroutine setup_sirius
    
end module sirius_interface