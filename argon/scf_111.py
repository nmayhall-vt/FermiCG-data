#!/usr/bin/env python

'''
Gamma point post-HF calculation needs only real integrals.
Methods implemented in finite-size system can be directly used here without
any modification.
'''

import numpy
from pyscf.pbc import gto, scf

#basis = "def2-SVPD"
#basis = "def2-SVP"
basis = "lanl2dz"

cell = gto.M(
    a = '''
    10.86000   0.00000   0.00000
    0.00000  10.86000   0.00000
    0.00000   0.00000  10.86000
    ''',
    atom = """
Ar   0.00000   0.00000   0.00000
Ar   0.00000   2.71500   2.71500
Ar   2.71500   0.00000   2.71500
Ar   2.71500   2.71500   0.00000
Ar   0.00000   0.00000   5.43000
Ar   0.00000   2.71500   8.14500
Ar   2.71500   0.00000   8.14500
Ar   2.71500   2.71500   5.43000
Ar   0.00000   5.43000   0.00000
Ar   0.00000   8.14500   2.71500
Ar   2.71500   5.43000   2.71500
Ar   2.71500   8.14500   0.00000
Ar   0.00000   5.43000   5.43000
Ar   0.00000   8.14500   8.14500
Ar   2.71500   5.43000   8.14500
Ar   2.71500   8.14500   5.43000
Ar   5.43000   0.00000   0.00000
Ar   5.43000   2.71500   2.71500
Ar   8.14500   0.00000   2.71500
Ar   8.14500   2.71500   0.00000
Ar   5.43000   0.00000   5.43000
Ar   5.43000   2.71500   8.14500
Ar   8.14500   0.00000   8.14500
Ar   8.14500   2.71500   5.43000
Ar   5.43000   5.43000   0.00000
Ar   5.43000   8.14500   2.71500
Ar   8.14500   5.43000   2.71500
Ar   8.14500   8.14500   0.00000
Ar   5.43000   5.43000   5.43000
Ar   5.43000   8.14500   8.14500
Ar   8.14500   5.43000   8.14500
Ar   8.14500   8.14500   5.43000
    """,
    basis = {'Ar': basis},
    ecp = {'Ar': basis},
    pseudo = {'Ar': basis}
    )


#kpts = cell.make_kpts([1,1,1])
#mf = scf.KRHF(cell, kpts = kpts).density_fit()

mf = scf.RHF(cell).density_fit()
mf.verbose = 4
mf.kernel()

#from pyscf.pbc.tools.k2gamma import k2gamma
#mf = k2gamma(mf)

#
# Import CC, TDDFT module from the molecular implementations
#
from pyscf import cc, tddft
#mycc = cc.CCSD(mf)
#mycc.kernel()

mytd = tddft.TDA(mf)
mytd.nstates = 5
mytd.verbose = 4
mytd.singlet = True 
mytd = mytd.run()
mytd.analyze()
