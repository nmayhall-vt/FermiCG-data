from functools import reduce
import numpy as np
import scipy

import pyscf
from pyscf import fci
from pyscf import gto, scf, ao2mo, lo, tdscf, cc


def tda_denisty_matrix(td, state_id):
    '''
    Taking the TDA amplitudes as the CIS coefficients, calculate the density
    matrix (in AO basis) of the excited states
    '''
    cis_t1 = td.xy[state_id][0]
    dm_oo =-np.einsum('ia,ka->ik', cis_t1.conj(), cis_t1)
    dm_vv = np.einsum('ia,ic->ac', cis_t1, cis_t1.conj())

    # The ground state density matrix in mo_basis
    mf = td._scf
    dm = np.diag(mf.mo_occ)

    # Add CIS contribution
    nocc = cis_t1.shape[0]
    # Note that dm_oo and dm_vv correspond to spin-up contribution. "*2" to
    # include the spin-down contribution
    dm[:nocc,:nocc] += dm_oo * 2
    dm[nocc:,nocc:] += dm_vv * 2

    # Transform density matrix to AO basis
    mo = mf.mo_coeff
    dm = np.einsum('pi,ij,qj->pq', mo, dm, mo.conj())
    return dm

mol.basis = '6-31g*'
mol.spin = 2
mol.charge = 0
mol.build()

#mf = scf.ROHF(mol).x2c()
mf = scf.UHF(mol)
mf.verbose = 4
mf.get_init_guess(mol, key='minao')
mf.conv_tol = 1e-9
#mf.level_shift = .1
#mf.diis_start_cycle = 4
#mf.diis_space = 10
mf.run(max_cycle=200)


n_triplets = 4
n_singlets = 4

avg_rdm1 = mf.make_rdm1()


# compute singlets
mytd = tdscf.TDA(mf)
mytd.singlet = True 
mytd = mytd.run(nstates=n_singlets)
mytd.analyze()
for i in range(mytd.nroots):
    avg_rdm1 += tda_denisty_matrix(mytd, i)

# compute triplets 
mytd = tdscf.TDA(mf)
mytd.singlet = False 
mytd = mytd.run(nstates=n_triplets)
mytd.analyze()
for i in range(mytd.nroots):
    avg_rdm1 += tda_denisty_matrix(mytd, i)

# normalize
avg_rdm1 = avg_rdm1 / (n_singlets + n_triplets + 1)


S = mf.get_ovlp()
np.save("overlap_mat", S)
np.save("density_mat", mf.make_rdm1())
np.save("rhf_mo_coeffs", mf.mo_coeff)
np.save("cis_sa_density_mat", avg_rdm1)


