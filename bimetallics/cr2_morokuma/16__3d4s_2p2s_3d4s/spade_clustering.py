import numpy as np


def spade_partitioning(O, U, Pv, S):
    """
    Find orbitals that most strongly overlap with the projector, P,  by doing O-O and V-V rotations. 
    [O,U] -> Of, Uf, Oe, Ue
    where Of (Uf) and Oe (Ue) are the Occupied (Unoccupied) orbitals of the fragment and environment, respectively.
    
    P[AO, frag]
    O[AO, occupied]
    U[AO, virtual]
    """

    print(" Partition %4i occupied and %4i virtual orbitals into a total of %4i orbitals" %(O.shape[1], U.shape[1], Pv.shape[1]))
    print(Pv.shape)
    PS = Pv.T @ S @ Pv

    P = Pv @ np.linalg.inv(PS) @ Pv.T

    #print(np.linalg.det(P.T @ S @ P))
    # assert(np.isclose(np.abs(np.linalg.det(P.T @ S @ P)), 1.0))
    nfrag = Pv.shape[1]

    _,so,Vo = np.linalg.svd(P @ S @ O, full_matrices=True)
    _,su,Vu = np.linalg.svd(P @ S @ U, full_matrices=True)

    s = np.concatenate((so, su))
    inds = [] 
    for i in range(len(so)):
        inds.append(i+1)
    for i in range(len(su)):
        inds.append(-i-1)
    inds = np.array(inds)

    perm = np.argsort(s)[::-1]
    inds = inds[perm]
    s = s[perm]

    print(" %16s %12s %12s" %("Singular Value", "Occupied", "Virtual"))
    for i in range(nfrag):
        print(" %16i %12.8f %12.8f" %(i, so[i], su[i]))



    print(" %16s %12s %12s" %("--", "--", "--"))
    print(" %16s %12s %12s" %("Index", "Sing. Val.", "Space"))
    for i in range(nfrag):
        ii = inds[i]
        if ii > 0:
            print(" %16i %12.8f %12s" %(i, s[i], "Occupied"))
        elif ii < 0:
            print(" %16i %12.8f %12s" %(i, s[i], "Virtual"))
        else:
            error("ArithmeticError")

    Of = O @ Vo.T
    Uf = U @ Vu.T

    
    return  Of, Uf