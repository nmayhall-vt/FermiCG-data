import pyscf
import pyscf.tools

from spade_clustering import *


molecule = """
H 0.0 0.0 0.0
H 2.0 0.0 0.0
H 4.0 0.0 0.0
He 0.0 2.0 0.0
He 2.0 2.0 0.0
He 4.0 2.0 0.0
He 6.0 2.0 0.0
He 8.0 2.0 0.0
"""

basis = "def2-svp"
pymol = pyscf.gto.Mole(
        atom    =   molecule,
        symmetry=   True,
        spin    =   1, # number of unpaired electrons
        charge  =   0,
        basis   =   basis)


pymol.build()
print("symmetry: ",pymol.topgroup)
# mf = pyscf.scf.UHF(pymol).x2c()
mf = pyscf.scf.ROHF(pymol)
mf.verbose = 4
mf.conv_tol = 1e-8
mf.conv_tol_grad = 1e-5
mf.chkfile = "scf.fchk"
mf.init_guess = "sad"

mf.run(max_cycle=200)

print(" Hartree-Fock Energy: %12.8f" % mf.e_tot)
# mf.analyze()

F = mf.get_fock()



# Find AO's corresponding to atoms 
full = []
frag1 = []
frag2 = []
frag3 = []
for ao_idx,ao in enumerate(mf.mol.ao_labels(fmt=False)):
    if ao[0] in (0, 1, 2, 3):
        if ao[2] in ("1s",):
            frag1.append(ao_idx)
            full.append(ao_idx)
    elif ao[0] in (4, 5):
        if ao[2] in ("1s", ):
            frag2.append(ao_idx)
            full.append(ao_idx)
    elif ao[0] == 6:
        if ao[2] in ("1s", ):
            frag3.append(ao_idx)
            full.append(ao_idx)


frags = [frag1, frag2, frag3]
print(frags)


C = mf.mo_coeff
S = mf.get_ovlp()
ndocc = mf.nelec[1]
nsing = mf.nelec[0] - ndocc
nvirt = mf.mol.nao - ndocc - nsing

# Just use alpha orbitals
Cdocc = mf.mo_coeff[:,0:ndocc]
Csing = mf.mo_coeff[:,ndocc:ndocc+nsing]
Cvirt = mf.mo_coeff[:,ndocc+nsing:ndocc+nsing+nvirt]

nbas = Cdocc.shape[0]

# Define projectors
Pfull = np.eye(nbas)[:,full] 

# Project MOs onto all fragments
# # Oact, Vact, Cenv, _ = spade_partitioning(Cdocc, Cvirt, Pfull, S)
# Cact = np.hstack((Oact,Vact))
[Oact, Sact, Vact], [Cenv, Cerr, _] = spade_partitioning((Cdocc, Csing, Cvirt), Pfull, S)
assert(Cerr.shape[1] == 0)
Cact = np.hstack((Oact,Vact))


# Project active orbitals onto fragments
init_fspace = []
Cfrags = []
for f in frags:
    Pf = np.eye(nbas)[:,f] 

    print()
    print(" Fragment: ", f)
    (Of, Sf, Vf), (_, _, _) = spade_partitioning((Cdocc, Csing, Cvirt), Pf, S)
    Cfrags.append(np.hstack((Of, Sf, Vf)))
    ndocc_f = Of.shape[1]
    init_fspace.append((ndocc_f+Sf.shape[1], ndocc_f))


# Orthogonalize Fragment orbitals
Cfrags = sym_ortho(Cfrags, S)

# Write Molden files for visualization
pyscf.tools.molden.from_mo(mf.mol, "Pfull.molden", Pfull)
pyscf.tools.molden.from_mo(mf.mol, "Cact.molden", Cact)
pyscf.tools.molden.from_mo(mf.mol, "Cenv.molden", Cenv)
for i in range(len(frags)):
    pyscf.tools.molden.from_mo(mf.mol, "Cfrag%i.molden"%i, Cfrags[i])

print(init_fspace)
