from functools import reduce
import numpy as np
import scipy

import pyscf
from pyscf import fci
from pyscf import gto, scf, ao2mo, lo, tdscf, cc


def tda_density_matrix(td, state_id):
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
V         -0.13019       -1.21502        0.60743
F          2.77625        0.38829        1.04287
F         -0.23723       -1.07820       -2.58664
F         -3.42500        3.39744        0.63807
F          3.59313       -3.42667       -1.02868
F         -1.86346        1.37527        1.44484
F          1.00229       -0.63544       -4.89564
F          2.77520        4.27680        3.66433
F          4.68942       -3.18108       -3.44751
F          3.80881        2.77571        1.66002
F         -5.74222        0.49019       -2.25012
F          0.58126        3.33516        5.04208
F          3.47665       -1.70425       -5.38711
F         -5.41584        2.97144       -1.17572
F         -0.50916        1.01017        4.41835
F         -4.05700       -1.49761       -1.61295
N         -2.06530       -1.17738        0.39300
N         -0.59933       -3.01506        1.67496
N          0.56294       -0.67517        2.35696
N          0.94156       -2.44915       -0.47758
N          0.38126        1.65740       -0.94694
C          2.17888        1.10783        2.01686
C          1.08077        0.56778        2.70484
C          1.61962       -2.25198       -1.67302
C         -3.58847        2.17949        0.09653
C          2.88957       -2.76692       -1.97592
C          2.72770        2.33406        2.32782
C          1.02675       -1.53838       -2.72412
C          1.64222       -1.33083       -3.93689
C          1.14186        2.60655        4.06809
C          0.12469       -2.95737        2.96634
H          1.05505       -3.27264        2.84321
H         -0.31806       -3.54867        3.62533
C         -2.77765        1.13449        0.48699
C         -2.73956       -2.46580        0.64868
H         -2.64218       -3.06719       -0.13165
H         -3.70431       -2.32262        0.81873
C          2.88037       -1.88195       -4.19690
C          2.22339        3.08864        3.36397
C          0.59190        1.38305        3.73038
C         -2.07229       -3.06993        1.86758
H         -2.32669       -2.56295        2.67884
H         -2.36386       -4.00905        1.98249
C          0.25223        0.64323       -0.41811
C         -2.89166       -0.15403       -0.04393
C         -4.59489        1.96734       -0.82396
C         -3.90790       -0.31280       -0.99491
C         -4.76004        0.71329       -1.36086
C         -0.13173       -4.15729        0.84936
H         -0.82147       -4.39166        0.17902
H          0.02193       -4.94746        1.42571
C          3.49383       -2.62249       -3.21126
C          0.13016       -1.52809        3.48023
H         -0.77568       -1.26690        3.78253
H          0.75560       -1.43808        4.24247
C          0.48758        2.93927       -1.65001
C          1.15632       -3.77325        0.14928
H          1.89974       -3.72684        0.80139
H          1.38468       -4.44569       -0.54049
C          1.15168        3.93514       -0.70262
H          0.61719        4.01370        0.11536
H          1.21153        4.81076       -1.13866
H          2.05129        3.62053       -0.47594
C          1.33404        2.70792       -2.89686
H          2.22191        2.38745       -2.63401
H          1.42343        3.54919       -3.39162
H          0.90035        2.03854       -3.46634
C         -0.92543        3.37805       -2.00909
H         -1.36421        2.67285       -2.52926
H         -0.88698        4.20032       -2.54091
H         -1.43477        3.54337       -1.18842
'''
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
dm1 = np.load("density_mat.npy")
mf.kernel(dm0=dm1)


nstates = 4


mytd = tdscf.TDA(mf)
mytd = mytd.run(nstates=nstates)
mytd.analyze()
for i in range(mytd.nroots):
    filename = "density_mat_cis_%i"%i
    np.save(filename, tda_density_matrix(mytd, i))



