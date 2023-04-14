using QCBase
using Printf
using ClusterMeanField
using FermiCG
using NPZ
using InCoreIntegrals
using RDM
using JLD2

pyscf = pyimport("pyscf");
fcidump = pyimport("pyscf.tools.fcidump");
ctx = fcidump.read("biphenyl-2mer-1d-fcidump");
h = ctx["H1"];
g = ctx["H2"];
ecore = ctx["ECORE"];
g = pyscf.ao2mo.restore("1", g, size(h,2))
ints = InCoreInts(ecore,h,g)
#
rdm1 = zeros(size(ints.h1))
#
#
na = 12
nb = 12
#
clusters_in    = [(1:6),(7:12),(13:18),(19:24)]
n_clusters = 4
#
## define clusters
cluster_list = [collect(1:6), collect(7:12), collect(13:18), collect(19:24)]
clusters = [MOCluster(i,collect(cluster_list[i])) for i = 1:length(cluster_list)]
init_fspace = [ (3,3) for i in 1:n_clusters]
display(clusters)

#run cmf_oo
e_cmf, U_cmf, d1  = ClusterMeanField.cmf_oo_diis(ints, clusters, init_fspace, RDM1(rdm1, rdm1), verbose=0, diis_start=3);

ints = orbital_rotation(ints,U_cmf)

M = 150

ref_fspace = FockConfig(init_fspace)
ecore = ints.h0

cluster_bases = FermiCG.compute_cluster_eigenbasis_spin(ints, clusters, d1, [3,3,3,3], ref_fspace, max_roots=M, verbose=1);
clustered_ham = FermiCG.extract_ClusteredTerms(ints, clusters)
cluster_ops = FermiCG.compute_cluster_ops(cluster_bases, ints);

FermiCG.add_cmf_operators!(cluster_ops, cluster_bases, ints, d1.a, d1.b);

nroots=17

ci_vector = FermiCG.TPSCIstate(clusters, ref_fspace, R=nroots)

ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig([1,1,1,1])] = zeros(Float64,nroots)

#Triplets
ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig([2,1,1,1])] = zeros(Float64,nroots)
ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig([1,2,1,1])] = zeros(Float64,nroots)
ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig([1,1,2,1])] = zeros(Float64,nroots)
ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig([1,1,1,2])] = zeros(Float64,nroots)

#Singlets
ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig([3,1,1,1])] = zeros(Float64,nroots)
ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig([1,3,1,1])] = zeros(Float64,nroots)
ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig([1,1,3,1])] = zeros(Float64,nroots)
ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig([1,1,1,3])] = zeros(Float64,nroots)

ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig([4,1,1,1])] = zeros(Float64,nroots)
ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig([1,4,1,1])] = zeros(Float64,nroots)
ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig([1,1,4,1])] = zeros(Float64,nroots)
ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig([1,1,1,4])] = zeros(Float64,nroots)

ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig([5,1,1,1])] = zeros(Float64,nroots)
ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig([1,5,1,1])] = zeros(Float64,nroots)
ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig([1,1,5,1])] = zeros(Float64,nroots)
ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig([1,1,1,5])] = zeros(Float64,nroots)

FermiCG.eye!(ci_vector)

ci_vector = FermiCG.add_spin_focksectors(ci_vector)

eci, v = FermiCG.tps_ci_direct(ci_vector, cluster_ops, clustered_ham);

e0a, v0a = FermiCG.tpsci_ci(v, cluster_ops, clustered_ham,
                            incremental  = true,
                            thresh_cipsi = 1e-3,
                            thresh_foi   = 1e-5,
                            thresh_asci  = -1, 
                            max_mem_ci = 100.0);

v0a = v0b
e0a = e0b

rotations = FermiCG.hosvd(v0a, cluster_ops)
for ci in clusters
    FermiCG.rotate!(cluster_ops[ci.idx], rotations[ci.idx])
    FermiCG.rotate!(cluster_bases[ci.idx], rotations[ci.idx])
    FermiCG.check_basis_orthogonality(cluster_bases[ci.idx])
end

tucker = [.8e-3, .6e-3]

for i in tucker
    e0b, v0b = FermiCG.tpsci_ci(ci_vector, cluster_ops, clustered_ham,
                                incremental  = true,
                                thresh_cipsi = i,
                                thresh_foi   = 1e-5,
                                thresh_asci  = -1,
                                max_mem_ci = 100.0);

    @time e2 = FermiCG.compute_pt2_energy(v0b, cluster_ops, clustered_ham, thresh_foi=1e-8);

    println()
    println("	*======TPSCI results======*")
    @printf("TCI Thresh: %8.6f  Dim:%8d\n",i,size(v0b)[1])
    println()
    @printf("TCI %5s %12s %12s\n", "Root", "E(0)", "E(2)") 
    for r in 1:nroots
        @printf("TCI %5s %12.8f %12.8f\n",r, e0b[r] + ecore, e0b[r] + e2[r] + ecore)
    end

    clustered_S2 = FermiCG.extract_S2(ci_vector.clusters)

    println()
    println("	*======TPSCI S2 results======*")
    @printf(" %-50s", "Compute FINAL S2 expectation values: ")
    @time s2 = FermiCG.compute_expectation_value_parallel(v0b, cluster_ops, clustered_S2)

    @printf(" %5s %12s %12s\n", "Root", "Energy", "S2") 
    for r in 1:nroots
        @printf(" %5s %12.8f %12.8f\n",r, e0b[r]+ecore, abs(s2[r]))
    end

    rotations = FermiCG.hosvd(v0b, cluster_ops)
    for ci in clusters
        FermiCG.rotate!(cluster_ops[ci.idx], rotations[ci.idx])
        FermiCG.rotate!(cluster_bases[ci.idx], rotations[ci.idx])
        FermiCG.check_basis_orthogonality(cluster_bases[ci.idx])
    end
end
