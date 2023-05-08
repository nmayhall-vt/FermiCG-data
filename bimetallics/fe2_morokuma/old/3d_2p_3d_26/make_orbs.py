import numpy as np
import pyscf
import pyscf.tools

molecule = """
 Fe 1.67785607 0.00052233 0.06475932
 Fe -1.67785607 -0.00052233 0.06475932
 O 0.00000000 0.00000000 -0.47099074
 Cl 1.87002704 -1.09796437 1.99091682
 Cl 2.93244917 -0.98210488 -1.47467288
 Cl 2.37160936 2.07954091 -0.50446591
 Cl -1.87002704 1.09796437 1.99091682
 Cl -2.93244917 0.98210488 -1.47467288
 Cl -2.37160936 -2.07954091 -0.50446591
"""

basis = "def2-svp"
pymol = pyscf.gto.Mole(
        atom    =   molecule,
        symmetry=   True,
        spin    =   10, # number of unpaired electrons
        charge  =   -2,
        basis   =   basis)


pymol.build()
print("symmetry: ",pymol.topgroup)
mf = pyscf.scf.ROHF(pymol)

Ccmf = np.load("Ccmf.npy")
pyscf.tools.molden.from_mo(mf.mol, "Ccmf.molden", Ccmf)
