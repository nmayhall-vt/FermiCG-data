using QCBase
using ClusterMeanField
using NPZ
using InCoreIntegrals
using RDM
using JLD2

C = npzread("mo_coeffs.npy")
h0 = npzread("ints_h0.npy")
h1 = npzread("ints_h1.npy")
h2 = npzread("ints_h2.npy")
ints = InCoreInts(h0, h1, h2)

clusters = [(1:6), (7:10), (11:16)]
init_fspace = [(3, 0), (4, 4), (3, 0)]

clusters = [(1,2,3,4,5,6,11,12,13,14,15,16), (7:10)]
init_fspace = [(6, 0), (4, 4)]

clusters = [(1:8), (9:16)]
init_fspace = [(5, 2), (2, 5)]

clusters = [(1,2,3,11,12,13),(4,5,6,14,15,16), (7:10)]
init_fspace = [(3, 3), (0,0), (4, 4)]

clusters = [(1, 2, 3),(4, ),(5, )]
init_fspace = [(3, 0), (1, 1), (1, 1)]# 
clusters = [(1, 2, 3),(4, 5)]
init_fspace = [(3, 0), (2, 2)]# 
clusters = [(1, 2, 3, 4, 5)]

init_fspace=  [(3, 2), (2, 2), (1, 1)]
clusters   =  [[1, 2, 3, 4], [5, 6], [7]]



# clusters = [(1:6)]
# init_fspace = [(6, 0)]


clusters = [MOCluster(i, collect(clusters[i])) for i = 1:length(clusters)]
display(clusters)

d1 = RDM1(n_orb(ints))


# # Do CMF
e_cmf, U, d1 = ClusterMeanField.cmf_oo(ints, clusters, init_fspace, RDM1(n_orb(ints)),
    verbose=0, sequential=true, max_iter_oo=50)

ints = orbital_rotation(ints, U)
C = C*U# Do CMF


e_cmf, U, d1 = ClusterMeanField.cmf_oo_diis(ints, clusters, init_fspace, d1,
                           maxiter_oo   = 500, 
                           maxiter_ci   = 200, 
                           maxiter_d1   = 200, 
                           verbose      = 0, 
                           tol_oo       = 1e-6, 
                           tol_d1       = 1e-9, 
                           tol_ci       = 1e-11, 
                           sequential   = true, 
                           alpha        = .1,
                           diis_start   = 1,
                           max_ss_size  = 24)

ints = orbital_rotation(ints, U)
C = C*U

npzwrite("Ccmf.npy", C)

@save "data_cmf.jld2" clusters init_fspace ints d1 e_cmf U 
