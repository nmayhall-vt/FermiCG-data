using QCBase
using ClusterMeanField
using FermiCG
using NPZ
using InCoreIntegrals
using RDM
using JLD2


@load "data_cmf.jld2"

M = 10 

init_space = [(4,4), (3,4), (4,4), (4,3)]

cluster_bases = FermiCG.compute_cluster_eigenbasis(ints, clusters, 
                                                   verbose=1, 
                                                   max_roots=M, 
                                                   init_fspace=init_fspace, 
                                                   rdm1a=d1.a, 
                                                   rdm1b=d1.b, 
                                                   T=Float64)

clustered_ham = FermiCG.extract_ClusteredTerms(ints, clusters)
cluster_ops = FermiCG.compute_cluster_ops(cluster_bases, ints);

FermiCG.add_cmf_operators!(cluster_ops, cluster_bases, ints, d1.a, d1.b);

nroots=4
ci_vector = FermiCG.TPSCIstate(clusters, ref_fock, R=nroots)


eci, v = tps_ci_direct(ci_vector, cluster_ops, clustered_ham);

ci_vector = add_spin_focksectors(ci_vector)
eci, v = tps_ci_direct(ci_vector, cluster_ops, clustered_ham);
