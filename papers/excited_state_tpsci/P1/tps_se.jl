using FermiCG, NPZ, JLD2
using PyCall
using LinearAlgebra
using Printf
using QCBase
using RDM
using ClusterMeanField


@load  "cmf_diis.jld2"

ref_fock = FockConfig(init_fspace)
ecore = ints.h0

M = 150

# Build Cluster basis
cluster_bases = FermiCG.compute_cluster_eigenbasis_spin(ints, clusters, d1, [3,3,3,3], ref_fock, max_roots=M, verbose=1);

# Build ClusteredOperator
clustered_ham = FermiCG.extract_ClusteredTerms(ints, clusters);

# Build Cluster Operators
cluster_ops = FermiCG.compute_cluster_ops(cluster_bases, ints);

# Add cmf hamiltonians for doing MP-style PT2 
FermiCG.add_cmf_operators!(cluster_ops, cluster_bases, ints, d1.a, d1.b, verbose=0);

nroots = 9
N = 4
Nk = 150

ci_vector = FermiCG.TPSCIstate(clusters, FermiCG.FockConfig(init_fspace), R=nroots);

ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig([1,1,1,1])] = zeros(Float64,nroots)
ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig([2,1,1,1])] = zeros(Float64,nroots)
ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig([1,2,1,1])] = zeros(Float64,nroots)
ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig([1,1,2,1])] = zeros(Float64,nroots)
ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig([1,1,1,2])] = zeros(Float64,nroots)
ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig([3,1,1,1])] = zeros(Float64,nroots)
ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig([1,3,1,1])] = zeros(Float64,nroots)
ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig([1,1,3,1])] = zeros(Float64,nroots)
ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig([1,1,1,3])] = zeros(Float64,nroots)

FermiCG.eye!(ci_vector)

e0, vec_out = FermiCG.tps_ci_davidson(ci_vector, cluster_ops, clustered_ham,
                        conv_thresh = 1e-4,
                        max_ss_vecs = 4)

display(e0)
display(ci_vector)
@time e2a = FermiCG.compute_pt2_energy(vec_out, cluster_ops, clustered_ham, thresh_foi=1e-8);
@printf("TCI %5s %12s %12s\n", "Root", "E(0)", "E(2)") 
for r in 1:nroots
    @printf("TCI %5s %12.8f %12.8f\n",r, e0[r] + ecore, e0[r] + e2a[r] + ecore)
end

@save  "tps_se_9roots.jld2" e0 vec_out e2a ecore

