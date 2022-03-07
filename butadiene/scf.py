from functools import reduce
import numpy as np
import scipy

import pyscf
from pyscf import fci
from pyscf import gto, scf, ao2mo, lo, tdscf, cc



def get_natural_orbital_active_space(rdm, S, thresh=.01):
   

    Ssqrt = scipy.linalg.sqrtm((S+S.T)/2.0)
    Sinvsqrt = scipy.linalg.inv(Ssqrt)

    print(" Number of electrons found %12.8f" %np.trace(S@rdm))

    Dtot = Ssqrt.T @ rdm @ Ssqrt
    #Dtot = Ssqrt.T @ ( da + db) @ Ssqrt
    D_evals, D_evecs = np.linalg.eigh((Dtot+Dtot.T)/2.0)

    sorted_list = np.argsort(D_evals)[::-1]
    D_evals = D_evals[sorted_list]
    D_evecs = D_evecs[:,sorted_list]

    act_list = []
    doc_list = []


    for idx,n in enumerate(D_evals):
        print(" %4i = %12.8f" %(idx,n),end="")
        if n < 2.0 - thresh:
            if n > thresh:
                act_list.append(idx)
                print(" Active")
            else:
                print(" Virt")
        else:
            doc_list.append(idx)
            print(" DOcc")

    print(" Number of active orbitals: ", len(act_list))
    print(" Number of doc    orbitals: ", len(doc_list))

    D_evecs = Sinvsqrt @ D_evecs
    Cdoc = D_evecs[:, doc_list]
    Cact = D_evecs[:, act_list]
    return Cdoc, Cact 




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

mol = gto.Mole()
mol.atom = '''
C         -2.30986        0.50909        0.01592
C         -0.98261        0.43259       -0.05975
C         -0.26676       -0.68753        0.57690
C          1.06323       -0.78274        0.51273
H         -2.84237        1.33139       -0.45081
H         -2.89518       -0.26464        0.20336
H         -0.43345        1.20285       -0.59462
H         -0.82592       -1.45252        1.10939
H          1.65570       -0.03776       -0.01031
H          1.57073       -1.61569        0.98847
'''


mol.basis = '6-31g*'
mol.spin = 0
mol.build()

#mf = scf.ROHF(mol).x2c()
mf = scf.RHF(mol)
mf.verbose = 4
mf.get_init_guess(mol, key='minao')
mf.conv_tol = 1e-9
#mf.level_shift = .1
#mf.diis_start_cycle = 4
#mf.diis_space = 10
mf.run(max_cycle=200)


#compute CCSD(T) to compare
mycc = pyscf.cc.CCSD(mf).run()
print('CCSD total energy', mycc.e_tot)
et = mycc.ccsd_t()
print('CCSD(T) total energy', mycc.e_tot + et)


n_triplets = 1
n_singlets = 1

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
print(" Number of electrons found %12.8f" %np.trace(S@avg_rdm1))

Cdoc, Cact = get_natural_orbital_active_space(avg_rdm1, S, thresh=.01)

# localize
Cact = pyscf.lo.PM(mol).kernel(Cact, verbose=4);
pyscf.tools.molden.from_mo(mol, "Cact.molden", Cact)

#
# Build integrals

# First get the density from the doubly occupied orbitals 
# to include in our effective 1 body operator
d1_embed = 2 * Cdoc @ Cdoc.T

h0 = pyscf.gto.mole.energy_nuc(mol)
h  = pyscf.scf.hf.get_hcore(mol)
j, k = pyscf.scf.hf.get_jk(mol, d1_embed, hermi=1)

h0 += np.trace(d1_embed @ ( h + .5*j - .25*k))

# Rotate 1electron terms to active space
h = Cact.T @ h @ Cact
j = Cact.T @ j @ Cact;
k = Cact.T @ k @ Cact;

h1 = h + j - .5*k;

# form 2e integrals in active space
nact = h.shape[0]
h2 = pyscf.ao2mo.kernel(mol, Cact, aosym="s4", compact=False)
h2.shape = (nact, nact, nact, nact)

np.save("integrals_h0", h0)
np.save("integrals_h1", h1)
np.save("integrals_h2", h2)
np.save("mo_coeffs_act", Cact)
np.save("mo_coeffs_doc", Cdoc)
np.save("overlap_mat", S)
np.save("density_mat", mf.make_rdm1())

