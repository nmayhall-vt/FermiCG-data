using QCBase
using ClusterMeanField
using FermiCG
using NPZ
using InCoreIntegrals
using RDM
using JLD2


@load "data_cmf.jld2"

M = 400

init_fspace = FockConfig(init_fspace)


cluster_bases = FermiCG.compute_cluster_eigenbasis_spin(ints, clusters, d1, [3,3], init_fspace, max_roots=M, verbose=1);

clustered_ham = FermiCG.extract_ClusteredTerms(ints, clusters)
cluster_ops = FermiCG.compute_cluster_ops(cluster_bases, ints);

FermiCG.add_cmf_operators!(cluster_ops, cluster_bases, ints, d1.a, d1.b);

nroots=6
ci_vector = FermiCG.TPSCIstate(clusters, init_fspace, R=nroots)

ci_vector = FermiCG.add_spin_focksectors(ci_vector)

eci, v = FermiCG.tps_ci_direct(ci_vector, cluster_ops, clustered_ham);

e0a, v0a = FermiCG.tpsci_ci(ci_vector, cluster_ops, clustered_ham, 
                            incremental  = true,
                            thresh_cipsi = 1e-4, 
                            thresh_foi   = 1e-6, 
                            thresh_asci  = -1);
    

rotations = FermiCG.hosvd(v0a, cluster_ops)
for ci in clusters
    FermiCG.rotate!(cluster_ops[ci.idx], rotations[ci.idx])
    FermiCG.rotate!(cluster_bases[ci.idx], rotations[ci.idx])
    FermiCG.check_basis_orthogonality(cluster_bases[ci.idx])
end


e0a, v0a = FermiCG.tpsci_ci(ci_vector, cluster_ops, clustered_ham, 
                            incremental  = true,
                            thresh_cipsi = 1e-4, 
                            thresh_foi   = 1e-6, 
                            thresh_asci  = -1);
    

