using FermiCG
using JLD2
using PyCall
#using Plots
using LinearAlgebra
using Printf
using QCBase
using RDM
using ClusterMeanField

@load  "cmf_diis.jld2"

M = 150

ref_fspace = FockConfig(init_fspace)
ecore = ints.h0

@load  "tucker_thresh_0.0004.jld2"

clustered_ham = FermiCG.extract_ClusteredTerms(ints, clusters)
cluster_ops = FermiCG.compute_cluster_ops(cluster_bases, ints);

FermiCG.add_cmf_operators!(cluster_ops, cluster_bases, ints, d1.a, d1.b);

nroots=17

FermiCG.clip!(v0b, thresh=0.0006)
e0, v11 = FermiCG.tps_ci_direct(v0b, cluster_ops, clustered_ham)
@time e2 = FermiCG.compute_pt2_energy(v11, cluster_ops, clustered_ham, thresh_foi=1e-8);

clustered_S2 = FermiCG.extract_S2(v11.clusters)
@time s2 = FermiCG.compute_expectation_value_parallel(v11, cluster_ops, clustered_S2)

@save "tucker_clip_0.0006.jld2" clusters cluster_bases e0 v11 e2 ecore s2

println()
println("	*======TPSCI results======*")
@printf("TCI Thresh: %8.6f  Dim:%8d\n",0.0006,size(v11)[1])
println()
@printf("TCI %5s %12s %12s\n", "Root", "E(0)", "E(2)") 
for r in 1:nroots
    @printf("TCI %5s %12.8f %12.8f\n",r, e0[r] + ecore, e0[r] + e2[r] + ecore)
end

println()
println("	*======TPSCI S2 results======*")
@printf(" %-50s", "Compute FINAL S2 expectation values: ")
@printf(" %5s %12s %12s\n", "Root", "Energy", "S2") 
for r in 1:nroots
    @printf(" %5s %12.8f %12.8f\n",r, e0[r]+ecore, abs(s2[r]))
end

FermiCG.clip!(v11, thresh=0.0008)
e0, v2 = FermiCG.tps_ci_direct(v11, cluster_ops, clustered_ham)
@time e2 = FermiCG.compute_pt2_energy(v2, cluster_ops, clustered_ham, thresh_foi=1e-8);

clustered_S2 = FermiCG.extract_S2(v2.clusters)
@time s2 = FermiCG.compute_expectation_value_parallel(v2, cluster_ops, clustered_S2)

@save "tucker_clip_0.0008.jld2" clusters cluster_bases e0 v2 e2 ecore s2

println()
println("	*======TPSCI results======*")
@printf("TCI Thresh: %8.6f  Dim:%8d\n",0.0008,size(v2)[1])
println()
@printf("TCI %5s %12s %12s\n", "Root", "E(0)", "E(2)") 
for r in 1:nroots
    @printf("TCI %5s %12.8f %12.8f\n",r, e0[r] + ecore, e0[r] + e2[r] + ecore)
end

println()
println("	*======TPSCI S2 results======*")
@printf(" %-50s", "Compute FINAL S2 expectation values: ")
@printf(" %5s %12s %12s\n", "Root", "Energy", "S2") 
for r in 1:nroots
    @printf(" %5s %12.8f %12.8f\n",r, e0[r]+ecore, abs(s2[r]))
end

