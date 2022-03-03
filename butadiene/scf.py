from functools import reduce
import numpy as np
import scipy

import pyscf
from pyscf import fci
from pyscf import gto, scf, ao2mo, lo, tdscf



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

Cdoc, Cact = get_natural_orbital_active_space(avg_rdm1, S)

# localize
Cact = pyscf.lo.PM(mol).kernel(Cact, verbose=4);
pyscf.tools.molden.from_mo(mol, "Cact.molden", Cact)


##
##   Build and store
#
#
#j, k = mf.get_jk()
#h    = mf.get_hcore()
#s = mf.get_ovlp()
#
#fa = h + j[0] + j[1] - k[0]
#fb = h + j[0] + j[1] - k[1]
#
#ca = mf.mo_coeff[0]
#cb = mf.mo_coeff[1]
#
#dm1 = mf.make_rdm1();
#da = dm1[0]
#db = dm1[1]
#print(np.max( (ca.T @ fa @ ca).diagonal() - mf.mo_energy[0]))
#print(np.max( (cb.T @ fb @ cb).diagonal() - mf.mo_energy[1]))
#
#print(np.max( (ca.T @ mf.get_fock()[0] @ ca).diagonal() - mf.mo_energy[0]))
#print(np.max( (cb.T @ mf.get_fock()[1] @ cb).diagonal() - mf.mo_energy[1]))
#
#np.save("uhf_mo_coeffs_a", mf.mo_coeff[0])
#np.save("uhf_mo_coeffs_b", mf.mo_coeff[1])
#np.save("uhf_j_a", j[0])
#np.save("uhf_j_b", j[1])
#np.save("uhf_k_a", k[0])
#np.save("uhf_k_b", k[1])
#np.save("uhf_h", h)
#np.save("uhf_f_a", mf.get_fock()[0])
#np.save("uhf_f_b", mf.get_fock()[1])
#np.save("uhf_rdm1_a", da)
#np.save("uhf_rdm1_b", db)
