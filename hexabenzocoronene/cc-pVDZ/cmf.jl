using QCBase
using ClusterMeanField
using NPZ
using InCoreIntegrals
using RDM
using JLD2
using Printf


h0 = npzread("/Users/arnab/arnab/workspace/cMF_data/new_pi_conj_sys/ints_h0.npy")
h1 = npzread("/Users/arnab/arnab/workspace/cMF_data/new_pi_conj_sys/ints_h1.npy")
h2 = npzread("/Users/arnab/arnab/workspace/cMF_data/new_pi_conj_sys/ints_h2.npy")
ints = InCoreInts(h0, h1, h2)

# Define clusters

clusters_in = [
    (1:6),   
    (7:12),  
    (13:18), 
    (19:24),
    (25:30),
    (31:36),
    (37:42)
]
clusters = [MOCluster(i,collect(clusters_in[i])) for i = 1:length(clusters_in)]

init_fspace = [
    (3,3),
    (3,3),
    (3,3),
    (3,3),
    (3,3),
    (3,3),
    (3,3)
];
print(clusters)
print(init_fspace)

rdm1 = RDM1(n_orb(ints))
# # # Do CMF
@time e_cmf, U, d1 = ClusterMeanField.cmf_oo(ints, clusters, init_fspace, rdm1,
    verbose=0, gconv=1e-8,tol_d1=1e-9,tol_ci=1e-11,sequential=true, max_iter_oo=100)
ints = orbital_rotation(ints, U)

@save "hexabenzocoronene.jld2" clusters init_fspace ints d1 e_cmf U
    