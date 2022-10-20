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
Pa = npzread("Pa.npy")
Pb = npzread("Pb.npy")
ints = InCoreInts(h0, h1, h2)
d1 = RDM1(Pa, Pb)

@printf(" Input energy:    %12.8f\n", compute_energy(ints, d1))
ssd1 = ssRDM1(d1)
ssd2 = ssRDM2(RDM2(d1))
@printf(" Input energy:    %12.8f\n", compute_energy(ints, ssd1, ssd2))
#clusters = [(1:4), (5:9), (10:14)]
#init_fspace = [(4,4), (5,0), (5,0)]
clusters = [(1:5), (6:10), (11:14)]
init_fspace = [(5,0), (0,5), (4,4)]


clusters = [MOCluster(i,collect(clusters[i])) for i = 1:length(clusters)]
display(clusters)

for ci in clusters
    rdmi = subset(d1, ci)
    @printf(" %4i alpha: %12.8f beta: %12.8f full: %12.8f\n", ci.idx, tr(rdmi.a), tr(rdmi.b), tr(d1))
end
d1 = RDM1(n_orb(ints))



#e_cmf, U, d1 = ClusterMeanField.cmf_oo_diis(ints, clusters, init_fspace, d1,
#                           maxiter_oo   = 200, 
#                           maxiter_ci   = 200, 
#                           maxiter_d1   = 200, 
#                           max_ss_size  = 8,
#                           verbose      = 0, 
#                           tol_oo       = 1e-7, 
#                           tol_d1       = 1e-9, 
#                           tol_ci       = 1e-11, 
#                           sequential   = false, 
#                           alpha        = .2,
#                           use_pyscf    = false,
#                           diis_start   = 300)
#
#ints = orbital_rotation(ints, U)

e_cmf, U, d1 = ClusterMeanField.cmf_oo_diis(ints, clusters, init_fspace, d1,
                           maxiter_oo   = 300, 
                           maxiter_ci   = 200, 
                           maxiter_d1   = 200, 
                           max_ss_size  = 12,
                           verbose      = 0, 
                           tol_oo       = 1e-7, 
                           tol_d1       = 1e-9, 
                           tol_ci       = 1e-11, 
                           sequential   = false, 
                           alpha        = .1,
                           use_pyscf    = false,
                           diis_start   = 10)

ints = orbital_rotation(ints, U)

C = npzread("mo_coeffs.npy")
S = npzread("overlap_mat.npy")
Ccmf = C*U
npzwrite("Ccmf.npy", Ccmf)

@save "data_cmf.jld2" clusters init_fspace ints d1 e_cmf U Ccmf
