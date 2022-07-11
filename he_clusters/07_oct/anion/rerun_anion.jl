using FermiCG
using PyCall
using Printf
using JLD2

@load "/home/nbraunsc/FermiCG-data/he_clusters/07_oct/anion_tpsci.jl.546735.scr/eq_after_cmf.jld2"

pyscf = pyimport("pyscf");
fcidump = pyimport("pyscf.tools.fcidump");
ctx = fcidump.read("/home/nbraunsc/FermiCG-data/he_clusters/07_oct/anion_tpsci.jl.546735.scr/fcidump.he07_oct");
h = ctx["H1"];
g = ctx["H2"];
ecore = ctx["ECORE"];

max_roots = 100

cluster_bases = FermiCG.compute_cluster_eigenbasis(ints, clusters, verbose=1, max_roots=max_roots, delta_elec=4, init_fspace=init_fspace, rdm1a=Da, rdm1b=Db);
#Build Clustered Operator
cluster_ham = FermiCG.extract_ClusteredTerms(ints, clusters);

#Build Cluster Operators
cluster_ops = FermiCG.compute_cluster_ops(cluster_bases, ints);
FermiCG.add_cmf_operators!(cluster_ops, cluster_bases, ints, Da, Db);

#Need to find reference state 
ref_fock = FermiCG.FockConfig(init_fspace)

cat1 = [(1,2),(1,1),(1,1),(1,1),(1,1),(1,1),(1,1)]
cat2 = [(1,1),(1,2),(1,1),(1,1),(1,1),(1,1),(1,1)]
cat3 = [(1,1),(1,1),(1,2),(1,1),(1,1),(1,1),(1,1)]
cat4 = [(1,1),(1,1),(1,1),(1,2),(1,1),(1,1),(1,1)]
cat5 = [(1,1),(1,1),(1,1),(1,1),(1,2),(1,1),(1,1)]
cat6 = [(1,1),(1,1),(1,1),(1,1),(1,1),(1,2),(1,1)]
cat7 = [(1,1),(1,1),(1,1),(1,1),(1,1),(1,1),(1,2)]

nroots = 7
#ci_vector = FermiCG.TPSCIstate(clusters, ref_fock, R=nroots)
ci_vector = FermiCG.TPSCIstate(clusters, FermiCG.FockConfig([(1,2),(1,1),(1,1),(1,1),(1,1),(1,1),(1,1)]), R=nroots);
#ci_vector = FermiCG.ClusteredState(clusters, ref_fock, R=nroots);
#Need to find the automated way to define these other excited configs away from ref state, example is to large
#to do by hand
#probably something to do with building p spaces and q spaces
CAT1 = FermiCG.FockConfig(cat1)
CAT2 = FermiCG.FockConfig(cat2)
CAT3 = FermiCG.FockConfig(cat3)
CAT4 = FermiCG.FockConfig(cat4)
CAT5 = FermiCG.FockConfig(cat5)
CAT6 = FermiCG.FockConfig(cat6)
CAT7 = FermiCG.FockConfig(cat7)

FermiCG.add_fockconfig!(ci_vector,CAT1)
FermiCG.add_fockconfig!(ci_vector,CAT2)
FermiCG.add_fockconfig!(ci_vector,CAT3)
FermiCG.add_fockconfig!(ci_vector,CAT4)
FermiCG.add_fockconfig!(ci_vector,CAT5)
FermiCG.add_fockconfig!(ci_vector,CAT6)
FermiCG.add_fockconfig!(ci_vector,CAT7)

ci_vector[CAT1][ClusterConfig([1,1,1,1,1,1,1])] = [1,0,0,0,0,0,0]
ci_vector[CAT2][ClusterConfig([1,1,1,1,1,1,1])] = [0,1,0,0,0,0,0]
ci_vector[CAT3][ClusterConfig([1,1,1,1,1,1,1])] = [0,0,1,0,0,0,0]
ci_vector[CAT4][ClusterConfig([1,1,1,1,1,1,1])] = [0,0,0,1,0,0,0]
ci_vector[CAT5][ClusterConfig([1,1,1,1,1,1,1])] = [0,0,0,0,1,0,0]
ci_vector[CAT6][ClusterConfig([1,1,1,1,1,1,1])] = [0,0,0,0,0,1,0]
ci_vector[CAT7][ClusterConfig([1,1,1,1,1,1,1])] = [0,0,0,0,0,0,1]

display(ci_vector)

thresh_list = [0.01, 0.001, 0.0001]

for thresh_cipsi in thresh_list
    e0, v0 = FermiCG.tpsci_ci(ci_vector, cluster_ops, cluster_ham,
                              #thresh_cipsi=1e-4, # Threshold for adding to P-space
                              thresh_cipsi=thresh_cipsi, # Threshold for adding to P-space
                              thresh_foi=1e-5,    # Threshold for keeping terms when defining FOIS
                              thresh_asci=1e-2,     # Threshold of P-space configs to search from
                              max_iter=10,
                              matvec=3);

    @time e2 = FermiCG.compute_pt2_energy(v0, cluster_ops, cluster_ham, thresh_foi=1e-8)
    name = "eq_tpsci_results"*string(thresh_cipsi)*".jld2"
    @save name e0 e2 v0 ecore

    #println()
    #println("	*======TPSCI results======*")
    #@printf("TCI Thresh: %8.6f  Dim:%8d\n",thresh_cipsi,size(v0)[1])
    #println()
    #@printf("TCI %5s %12s %12s\n", "Root", "E(0)", "E(2)") 
    #for r in 1:nroots
    #    @printf("TCI %5s %12.8f %12.8f\n",r, e0[r] + ecore, e0[r] + e2[r] + ecore)
    #end

    for r in 1:nroots
        @printf("TPSCI %5s %12.8f %12.8f\n",r, e0[r] + ecore, e0[r] + e2[r] + ecore)
        display(v0,thresh=1e-4,root=r)
    end

    thresh_print = 1e-2
    for r in 1:nroots
        println("\n ### Root: ", r, " ###")
        idx = 1
        for (fock,configs) in v0.data
            length(v0.clusters) == length(fock) || throw(Exception)
            length(v0.data[fock]) > 0 || continue
            @printf(" Dim %4i fock_space: ",length(v0.data[fock]))
            [@printf(" %-2i(%i:%i) ",fii,fi[1],fi[2]) for (fii,fi) in enumerate(fock)]
            println()

            for (config, value) in v0.data[fock]
                if abs(value[r]) > thresh_print
                    @printf(" %5i",idx)
                    for c in config
                        @printf("%3i",c)
                    end
                    @printf(":%12.5f\n",value[r])
                end
                idx += 1
            end
        end
    end

    println()
    println("	*======TPSCI results======*")
    #@printf("TCI Thresh: %8.6f  Dim:%8d\n",1e-3,size(v0)[1])
    @printf("TCI Thresh: %8.6f  Dim:%8d\n",0.001,size(v0)[1])
    println()
    @printf("TCI %5s %12s %12s\n", "Root", "E(0)", "E(2)") 
    for r in 1:nroots
        @printf("TCI %5s %12.8f %12.8f\n",r, e0[r] + ecore, e0[r] + e2[r] + ecore)
    end
end





