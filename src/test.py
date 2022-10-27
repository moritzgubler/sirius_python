#!/usr/bin/env python3
# from symbol import parameters


def setup_sirius(**kwargs):

    required_param_set = {
        "gk_cutoff",
        'ngridk',
        'pw_cutoff',
        'shift_k',
        "xc_functionals"
    }
    
    param_list = {"auto_rmt",
        "aw_cutoff",
        "core_relativity",
        "density_tol",
        "electronic_structure_method", 
        "energy_tol",
        "extra_charge",
        "gamma_point",
        "gk_cutoff",
        "hubbard_correction",
        "lmax_apw",
        "lmax_pot",
        "lmax_rho",
        "molecule",
        "ngridk",
        "nn_radius",
        "num_bands",
        "num_dft_iter",
        "num_fv_states",
        "num_mag_dims",
        "precision_gs",
        "precision_hs",
        "precision_wf",
        "pw_cutoff",
        "reduce_aux_bf",
        "shiftk",
        "smearing",
        "smearing_width",
        "so_correction",
        "use_ibz",
        "use_scf_correction",
        "use_symmetry",
        "valence_relativity",
        "vk",
        "xc_dens_tre",
        "xc_functionals"
    }

    for par in required_param_set:
        if par not in kwargs:
            print(par, 'not given')

    json_dict = {'parameters': {}}
    for kw in kwargs:
        if kw in param_list:
            print(kw)
            json_dict['parameters'][kw] = kwargs[kw]
            

    print(json_dict)


if __name__ == "__main__":
    setup_sirius(gk_cutoff = 3.0, pw_cutoff=5.0)