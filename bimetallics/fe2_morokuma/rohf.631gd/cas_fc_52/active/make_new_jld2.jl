using QCBase
using ClusterMeanField
using NPZ
using InCoreIntegrals
using RDM
using JLD2
using Printf


@load "../data_cmf.jld2"

C = npzread("mo_coeffs.npy")
h0 = npzread("ints_h0.npy")
h1 = npzread("ints_h1.npy")
h2 = npzread("ints_h2.npy")

ints = InCoreInts(h0, h1, h2)


init_fspace =  [(5, 0), (3, 3), (5, 0), (3, 3), (3, 3), (3, 3), (3, 3), (3, 3), (3, 3)]
clusters    =  [[16, 17, 18, 19, 20], [21, 22, 23, 24, 25, 26], [27, 28, 29, 30, 31], [32, 33, 34, 35, 36, 37], [38, 39, 40, 41, 42, 43], [44, 45, 46, 47, 48, 49], [50, 51, 52, 53, 54, 55], [56, 57, 58, 59, 60, 61], [62, 63, 64, 65, 66, 67]]


active_indices = []


for i in 1:length(clusters)
    for j in 1:length(clusters[i])
        push!(active_indices, clusters[i][j])
        clusters[i][j] -= 15
    end
end


d1 = RDM1(d1.a[active_indices, active_indices], d1.b[active_indices, active_indices])

clusters = [MOCluster(i, collect(clusters[i])) for i = 1:length(clusters)]

@save "data_cmf.jld2" clusters init_fspace ints d1 e_cmf U 
