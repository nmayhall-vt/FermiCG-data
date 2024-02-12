import numpy as np
import time
from math import factorial
from pyscf_helper import *
import pyscf
from pyscf import gto, scf, ao2mo,  lo
from pyscf.tools import molden
ttt = time.time()
def reorder_integrals(idx,h,g):
# {{{
    h = h[:,idx] 
    h = h[idx,:] 

    g = g[:,:,:,idx] 
    g = g[:,:,idx,:] 
    g = g[:,idx,:,:] 
    g = g[idx,:,:,:] 
    return h,g
# }}}

###     PYSCF INPUT
r0 = 2.0
molecule = '''
C 4.93354 -2.84833 0.00000
C 4.93321 -1.45861 0.00000
C 3.73075 -0.73168 0.00000
C 2.49032 -1.43783 0.00000
C 1.23513 -0.71311 0.00000
C 1.23513 0.71311 0.00000
C 2.49032 1.43783 0.00000
C 3.73075 0.73168 0.00000
C 4.93321 1.45861 0.00000
C 4.93354 2.84833 0.00000
C 3.72991 3.54302 0.00000
C 2.49908 2.86514 0.00000
C 1.23174 3.59686 -0.00000
C 1.20340 5.00170 -0.00000
C 0.00000 5.69682 -0.00000
C -1.20340 5.00170 -0.00000
C -1.23174 3.59686 -0.00000
C 0.00000 2.87564 -0.00000
C 0.00000 1.42626 -0.00000
C -1.23513 0.71311 -0.00000
C -2.49032 1.43783 0.00000
C -2.49908 2.86514 -0.00000
C -3.72991 3.54302 -0.00000
C -4.93354 2.84833 0.00000
C -4.93321 1.45861 0.00000
C -3.73075 0.73168 0.00000
C -3.73075 -0.73168 0.00000
C -4.93321 -1.45861 0.00000
C -4.93354 -2.84833 0.00000
C -3.72991 -3.54302 0.00000
C -2.49908 -2.86514 -0.00000
C -2.49032 -1.43783 -0.00000
C -1.23513 -0.71311 -0.00000
C -0.00000 -1.42626 -0.00000
C -0.00000 -2.87564 -0.00000
C -1.23174 -3.59686 -0.00000
C -1.20340 -5.00170 -0.00000
C 0.00000 -5.69682 -0.00000
C 1.20340 -5.00170 -0.00000
C 1.23174 -3.59686 -0.00000
C 2.49908 -2.86514 -0.00000
C 3.72991 -3.54302 0.00000
H 5.87949 -3.39439 0.00000
H 5.89025 -0.93984 0.00000
H 5.89025 0.93984 0.00000
H 5.87949 3.39439 0.00000
H 3.75925 4.63123 0.00000
H 2.13117 5.57119 -0.00000
H 0.00000 6.78906 -0.00000
H -2.13117 5.57119 -0.00000
H -3.75925 4.63123 -0.00000
H -5.87949 3.39439 0.00000
H -5.89025 0.93984 0.00000
H -5.89025 -0.93984 0.00000
H -5.87949 -3.39439 0.00000
H -3.75925 -4.63123 0.00000
H -2.13117 -5.57119 -0.00000
H 0.00000 -6.78906 -0.00000
H 2.13117 -5.57119 -0.00000
H 3.75925 -4.63123 -0.00000
'''
charge = 0
spin  = 0
basis_set = 'cc-pVDZ'

npoly = 21
na = npoly
nb = npoly

orb_basis = 'PM'
cas_nel = 2*npoly
cas_norb = 2*npoly
pyscf.lib.num_threads(1)  #with degenerate states and multiple processors there can be issues
#PYSCF inputs

mol = gto.Mole()
mol.atom = molecule

mol.max_memory = 1000 # MB
mol.symmetry = True
mol.charge = charge
mol.spin = spin
mol.basis = basis_set
mol.build()
print("symmetry")
print(mol.topgroup)

#SCF 

#mf = scf.RHF(mol).run(init_guess='atom')
mf = scf.RHF(mol)
mf.verbose = 4
mf.conv_tol = 1e-8
mf.conv_tol_grad = 1e-5
mf.chkfile = "scf.fchk"
mf.init_guess = "minao"
mf.run()
#C = mf.mo_coeff #MO coeffs
enu = mf.energy_nuc()


h,ecore,g,C = get_pi_space(mol,mf,cas_norb,cas_nel,local=True)

idx = e1_order(h,1e-1) 
h,g = reorder_integrals(idx,h,g)
print(h)
C = C[:,idx] # make sure u reorder this too
clusters = [ [0,1,2,5,6,11], [3,4,8,10,13,16],[7,9,12,15,18,21],[14,19,20,25,26,31],[22,27,28,32,33,36],[24,29,30,34,17,23],[35,37,38,39,40,41]]
init_fspace = [(3,3),(3,3),(3,3),(3,3),(3,3),(3,3),(3,3)]
block= [0,1,2,5,6,11, 3,4,8,10,13,16,7,9,12,15,18,21,14,19,20,25,26,31,22,27,28,32,33,36,24,29,30,34,17,23,35,37,38,39,40,41]
h,g = reorder_integrals(block,h,g)
C = C[:,block]

np.save("ints_h0_pi_cg_4", ecore)
np.save("ints_h1_pi_cg_4", h)
np.save("ints_h2_pi_cg_4", g)
np.save("mo_coeffs_pi_cg_4", C)
    
molden.from_mo(mol, 'cas_hbc.molden', C)
    
