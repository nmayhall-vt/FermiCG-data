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

clusters = [(1:8), (9:18), (19:26), (27:36)]
init_fspace = [(4,4), (2,5), (4,4), (5,2)]
#init_fspace = [(1,1), (2,0), (1,1), (2,0)]

clusters = [MOCluster(i,collect(clusters[i])) for i = 1:length(clusters)]
display(clusters)

d1 = RDM1(n_orb(ints))

ecmf, d1dict, d2dict = cmf_ci(ints, clusters, init_fspace, d1,
                    verbose=2, sequential=true, dconv=1e-5)


e_cmf, U, d1 = cmf_oo(ints, clusters, init_fspace, d1,
                           verbose=0, gconv=1e-5, method="bfgs",sequential=true)

ints = orbital_rotation(ints, U)

@save "data_cmf.jld2" clusters init_fspace ints d1 e_cmf U 
