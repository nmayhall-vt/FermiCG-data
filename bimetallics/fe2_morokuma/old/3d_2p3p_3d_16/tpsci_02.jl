
rotations = FermiCG.hosvd(v0b, cluster_ops)
for ci in clusters
    FermiCG.rotate!(cluster_ops[ci.idx], rotations[ci.idx])
    FermiCG.rotate!(cluster_bases[ci.idx], rotations[ci.idx])
    FermiCG.check_basis_orthogonality(cluster_bases[ci.idx])
end


e0b, v0b = FermiCG.tpsci_ci(ci_vector, cluster_ops, clustered_ham, 
                            incremental  = true,
                            thresh_cipsi = .2e-3, 
                            thresh_foi   = 1e-6, 
                            thresh_asci  = -1);
    

