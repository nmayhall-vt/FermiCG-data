using QCBase
using ClusterMeanField 
using NPZ
using InCoreIntegrals
using RDM
using JLD2

C = npzread("mo_coeffs.npy")
S = npzread("overlap_mat.npy")
h0 = npzread("ints_h0.npy")
h1 = npzread("ints_h1.npy")
h2 = npzread("ints_h2.npy")
ints = InCoreInts(h0, h1, h2)

clusters = [(1:12), (13:20), (21:32)]
init_fspace = [(6,1), (4,4), (1,6)]

clusters = [MOCluster(i,collect(clusters[i])) for i = 1:length(clusters)]
display(clusters)

d1 = RDM1(n_orb(ints))


e_cmf, U, d1 = ClusterMeanField.cmf_oo_diis(ints, clusters, init_fspace, d1,
                           maxiter_oo   = 200, 
                           maxiter_ci   = 200, 
                           maxiter_d1   = 200, 
                           verbose      = 0, 
                           tol_oo       = 1e-7, 
                           tol_d1       = 1e-9, 
                           tol_ci       = 1e-11, 
                           sequential   = true, 
                           alpha        = .1,
                           diis_start   = 1)

ints = orbital_rotation(ints, U)

C = npzread("mo_coeffs.npy")
S = npzread("overlap_mat.npy")
Ccmf = C*U
npzwrite("Ccmf.npy", Ccmf)

@save "data_cmf.jld2" clusters init_fspace ints d1 e_cmf U Ccmf
