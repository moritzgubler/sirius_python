import sirius
from ase.io import read
import numpy as np
angbohr = 1.8897161646320724
at = read('si.ascii')
pos = angbohr * at.get_positions().T
lat = angbohr * at.get_cell(True).T
nat = at.get_global_number_of_atoms()
symb = at.get_chemical_symbols()


at_names = np.empty((nat, 2), dtype='c')
for i in range(nat):
    at_names[i] = str(symb[i])

print(np.shape(pos))

kgrid = 3*np.ones(3, dtype=np.int32).T
kshift = np.zeros(3, dtype=np.int32).T

at_names = np.array(symb, dtype='c').T

json_dir = './'
pw_cutoff = 500.0
gk_cutoff = 90.0

print(type(pos))


sirius.mod_sirius_energy.setup_sirius(pos, lat, at_names, json_dir, pw_cutoff, gk_cutoff, kgrid, kshift, nat)

fxyz,etot,deralat = sirius.mod_sirius_energy.energyandforces_sirius(pos, lat, nat)

sirius.mod_sirius_energy.exit_sirius()
sirius.mod_sirius_energy.shutdown_sirius_interface()

print(etot)




# print('hello from python')
# 
# test.hello_world_mpi()
# 
# print('I should get 8 processes')
# 
# test.test_print()
# 
# print('now only one process')
# 
# test.release_process()
# 
# print('asdf')
# 
# test.finalize()