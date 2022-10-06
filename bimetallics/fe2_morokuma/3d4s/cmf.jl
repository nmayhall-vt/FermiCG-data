using QCBase
using ClusterMeanField 
using NPZ
using InCoreIntegrals
using RDM
using JLD2

h0 = npzread("ints_h0.npy")
h1 = npzread("ints_h1.npy")
h2 = npzread("ints_h2.npy")
ints = InCoreInts(h0, h1, h2)

clusters = [(1:12), (13:24)]
init_fspace = [(6,1), (1,6)]

clusters = [MOCluster(i,collect(clusters[i])) for i = 1:length(clusters)]
display(clusters)

d1 = RDM1(n_orb(ints))


e_cmf, U, d1 = ClusterMeanField.cmf_oo_diis(ints, clusters, init_fspace, d1,
                           max_iter_oo=200, verbose=0, conv_oo=1e-6, conv_ci=1e-9, sequential=true, alpha=.1)

ints = orbital_rotation(ints, U)

@save "data_cmf.jld2" clusters init_fspace ints d1 e_cmf U 
