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





mol = gto.Mole()
mol.atom = '''
H           -3.426100        -2.240400         5.488400
H           -5.627400        -1.077000         5.214700
C           -3.653500        -1.732700         4.551600
H           -1.767100        -2.237000         3.663900
C           -4.907300        -1.068800         4.394700
H           -6.163100         0.096400         3.101400
C           -2.725800        -1.732100         3.540600
H           -0.300300         1.083200        -5.235700
C           -5.209800        -0.419000         3.224900
C           -2.996100        -1.063600         2.307300
H           -1.103000        -1.532900         1.397700
H           -0.427000        -0.802900        -0.856600
H            0.236100        -0.097900        -3.127300
C           -1.019300         1.073000        -4.415000
H           -2.498800         2.251900        -5.503400
C           -4.274000        -0.392400         2.144500
H           -5.501500         0.794400         0.831000
C           -2.061300        -1.027200         1.271800
C           -1.382000        -0.289500        -0.977200
C           -0.717100         0.418000        -3.247600
C           -2.272000         1.739500        -4.569000
H           -4.157600         2.241200        -3.678700
C           -4.546300         0.281700         0.953400
C           -2.324300        -0.340200         0.070400
C           -1.652800         0.387400        -2.167000
C           -3.199800         1.734100        -3.558400
C           -3.604400         0.330900        -0.094300
C           -2.930200         1.059100        -2.329200
C           -3.866500         1.018700        -1.295500
H           -4.824300         1.525600        -1.421700
H            6.954300         5.020900        -5.781400
H            9.130800         5.788500        -4.794600
C            7.130300         4.913900        -4.709900
H            5.199700         4.047800        -4.322900
C            8.368100         5.351000        -4.148200
H            9.543000         5.570100        -2.368400
C            6.155900         4.370100        -3.908500
H            4.469100         3.270900        -2.103400
C            8.601900         5.233600        -2.800800
C            6.366300         4.207300        -2.503300
C            5.408400         3.620800        -1.673100
C            7.624000         4.655100        -1.933900
H            8.793800         4.849500        -0.136000
H            3.726000         2.535700         0.128600
C            5.630800         3.465600        -0.289200
C            7.849900         4.509800        -0.563800
C            4.668600         2.877900         0.557100
C            6.887300         3.923700         0.283000
H            8.046600         4.121200         2.098700
H            2.974600         1.818600         2.362500
C            4.893400         2.733400         1.927600
C            7.108600         3.769300         1.667200
C            3.915400         2.155900         2.795000
H            3.383500         1.608200         4.789600
C            6.150200         3.183400         2.497000
H            7.316000         3.346600         4.316500
C            4.148200         2.042400         4.143100
C            6.360000         3.023000         3.902600
C            5.385400         2.480900         4.704600
H            5.560600         2.377300         5.776700
H            1.028200         5.020900        -5.781400
H            3.204700         5.788500        -4.794600
C            1.204300         4.913900        -4.709900
H           -0.726400         4.047800        -4.322900
C            2.442000         5.351000        -4.148200
H            3.616900         5.570100        -2.368400
C            0.229800         4.370100        -3.908500
H           -1.457000         3.270900        -2.103400
C            2.675800         5.233600        -2.800800
C            0.440300         4.207300        -2.503300
C           -0.517600         3.620800        -1.673100
C            1.697900         4.655100        -1.933900
H            2.867800         4.849500        -0.136000
H           -2.200000         2.535700         0.128600
H           -2.951500         1.818600         2.362500
C           -0.295200         3.465600        -0.289200
C            1.923900         4.509800        -0.563800
H           -2.542500         1.608200         4.789600
C           -1.257500         2.877900         0.557100
C           -2.010600         2.155900         2.795000
C            0.961300         3.923700         0.283000
H            2.120600         4.121200         2.098700
C           -1.777900         2.042400         4.143100
H           -0.365500         2.377300         5.776700
C           -1.032700         2.733400         1.927600
C            1.182500         3.769300         1.667200
C           -0.540600         2.480900         4.704600
C            0.224100         3.183400         2.497000
C            0.434000         3.023000         3.902600
H            1.389900         3.346600         4.316500
H            3.427300         2.251900        -5.503400
H            5.625800         1.083200        -5.235700
C            3.654100         1.739500        -4.569000
H            1.768500         2.241200        -3.678700
H            0.298700        -1.077000         5.214700
C            4.906700         1.073000        -4.415000
H            6.162100        -0.097900        -3.127300
C            2.726300         1.734100        -3.558400
H           -0.237000         0.096400         3.101400
C            1.018700        -1.068800         4.394700
H            2.500000        -2.240400         5.488400
H            0.424600         0.794400         0.831000
H            1.101800         1.525600        -1.421700
C            5.209000         0.418000        -3.247600
C            2.995900         1.059100        -2.329200
C            0.716200        -0.419000         3.224900
C            2.272600        -1.732700         4.551600
H            4.159000        -2.237000         3.663900
C            1.379700         0.281700         0.953400
C            2.059600         1.018700        -1.295500
C            4.273300         0.387400        -2.167000
H            5.499100        -0.802900        -0.856600
C            1.652100        -0.392400         2.144500
C            3.200300        -1.732100         3.540600
C            2.321700         0.330900        -0.094300
C            4.544100        -0.289500        -0.977200
C            2.930000        -1.063600         2.307300
C            3.601700        -0.340200         0.070400
C            3.864800        -1.027200         1.271800
H            4.823000        -1.532900         1.397700
'''

np.save("xyz.npy", mol.atom)

mol.basis = '6-31g*'
mol.spin = 0
mol.build()

mf = scf.RHF(mol).density_fit()
mf.verbose = 4
mf.get_init_guess(mol, key='minao')
mf.conv_tol = 1e-9

# load precomputed data
C = np.load("../rhf_mo_coeffs.npy")
avg_rdm1 = np.load("../cis_sa_density_mat.npy")

S = mf.get_ovlp()

print(avg_rdm1.shape)
print(S.shape)

print(" Number of electrons found %12.8f" %np.trace(S@avg_rdm1))

Cdoc, Cact = get_natural_orbital_active_space(avg_rdm1, S, thresh=.00275)

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


