using FermiCG
using PyCall
using Plots
using LinearAlgebra
using Printf
using JLD2

molecule = "
He 0.00000000 0.00000000 0.00000000
He 2.82842713 0.00000000 0.00000000
He 0.00000000 2.82842713 0.00000000
He 2.82842712 2.82842712 0.00000000
He 1.41421356 1.41421356 2.00000000
He 1.41421356 1.41421356 -2.00000000
He 1.41421356 1.41421356 0.00000000
"
atoms = []
for (li,line) in enumerate(split(rstrip(lstrip(molecule)), "\n"))
    l = split(line)
    push!(atoms, Atom(li, l[1], parse.(Float64,l[2:4])))
end

#basis = "aug-cc-pvdz" #9 orbs on each He
basis = "cc-pvdz" #5 orbs on each He
#basis = "cc-pvtz" # 14 orbs on each He
#basis = "sto-3g"

# Create FermiCG.Molecule type
mol = Molecule(0,1,atoms,basis)

pyscf = pyimport("pyscf")
pymol = pyscf.gto.Mole(atom=molecule, spin=0, charge=0, basis=basis)
pymol.build()
mf = pyscf.scf.RHF(pymol).run()
s = mf.get_ovlp(pymol)

lo = pyimport("pyscf.lo.orth")
lo_ao = lo.lowdin(s)
println("size of Lowdin ortho AO's:", size(lo_ao))

FermiCG.pyscf_write_molden(mol, lo_ao, filename="lowdin_ao_ccpvdz.molden")

#write fci dump file from the modified mo coefficients
tools = pyimport("pyscf.tools")
tools.fcidump.from_mo(pymol, "fcidump.he07_oct", lo_ao)

#Can just read in pyscf dump file for integrals (once you have already run an scf calculation)
pyscf = pyimport("pyscf");
fcidump = pyimport("pyscf.tools.fcidump");
ctx = fcidump.read("fcidump.he07_oct");
h = ctx["H1"];
g = ctx["H2"];
ecore = ctx["ECORE"];
g = pyscf.ao2mo.restore("1", g, size(h,2))

#This one below was not working. Error: setfield! immutable struct of type InCoreInts cannot be changed
ints = InCoreInts(ecore,h,g);


#Run cmf
#Define clusters and intial Fock space for inital CMF calc for 14 orbs each He
#clusters_in = [(1:14),(15:28), (29:42), (43:56), (57:70), (71:84), (85:98)]

#Define clusters and intial Fock space for inital CMF calc for 9 orbs each He
#clusters_in = [(1:9),(10:18), (19:27), (28:36), (37:45), (46:54), (55:63)]

#Define clusters and intial Fock space for inital CMF calc for 5 orbs each He
clusters_in = [(1:5),(6:10), (11:15), (16:20), (21:25), (26:30), (31:35)]
init_fspace = [(1,1),(1,1),(1,1),(1,1),(1,1),(1,1),(1,1)]
rdm1 = zeros(size(ints.h1))
na=7
nb=7

#Define clusters now using FermiCG code
clusters = [Cluster(i,collect(clusters_in[i])) for i = 1:length(clusters_in)]
display(clusters)

@save "before_cmf.jld2" ints clusters init_fspace

print(size(ints.h1))
rdm1 = zeros(size(ints.h1))

#do a CMF calculation to optimize cluster orbitals

e_cmf, U, Da, Db = FermiCG.cmf_oo(ints, clusters, init_fspace, rdm1, rdm1, max_iter_oo=200, verbose=0, gconv=1e-6, method="bfgs");


#rotate the integrals by the cmf calculation
ints = FermiCG.orbital_rotation(ints, U);
max_roots = 100

#Build Cluster Basis (delta n is here)
cluster_bases = FermiCG.compute_cluster_eigenbasis(ints, clusters, verbose=1, max_roots=max_roots, init_fspace=init_fspace, rdm1a=Da, rdm1b=Db);

@save "eq_after_cmf.jld2" ints Da Db e_cmf cluster_bases clusters init_fspace

FermiCG.pyscf_write_molden(mol,lo_ao*U, filename="cmf_he07_oct.molden");

#Build Clustered Operator
cluster_ham = FermiCG.extract_ClusteredTerms(ints, clusters);

#Build Cluster Operators
cluster_ops = FermiCG.compute_cluster_ops(cluster_bases, ints);
FermiCG.add_cmf_operators!(cluster_ops, cluster_bases, ints, Da, Db);

cat1 = [(1,0),(1,1),(1,1),(1,1),(1,1),(1,1),(1,1)]
cat2 = [(1,1),(1,0),(1,1),(1,1),(1,1),(1,1),(1,1)]
cat3 = [(1,1),(1,1),(1,0),(1,1),(1,1),(1,1),(1,1)]
cat4 = [(1,1),(1,1),(1,1),(1,0),(1,1),(1,1),(1,1)]
cat5 = [(1,1),(1,1),(1,1),(1,1),(1,0),(1,1),(1,1)]
cat6 = [(1,1),(1,1),(1,1),(1,1),(1,1),(1,0),(1,1)]
cat7 = [(1,1),(1,1),(1,1),(1,1),(1,1),(1,1),(1,0)]


#Need to find reference state 
ref_fock = FermiCG.FockConfig(init_fspace)
CAT1 = FermiCG.FockConfig(cat1)
CAT2 = FermiCG.FockConfig(cat2)
CAT3 = FermiCG.FockConfig(cat3)
CAT4 = FermiCG.FockConfig(cat4)
CAT5 = FermiCG.FockConfig(cat5)
CAT6 = FermiCG.FockConfig(cat6)
CAT7 = FermiCG.FockConfig(cat7)

nroots = 7
#ci_vector = FermiCG.TPSCIstate(clusters, ref_fock, R=nroots)
ci_vector = FermiCG.TPSCIstate(clusters, FermiCG.FockConfig([(1,0),(1,1),(1,1),(1,1),(1,1),(1,1),(1,1)]), R=nroots);
#ci_vector = FermiCG.ClusteredState(clusters, ref_fock, R=nroots);
#Need to find the automated way to define these other excited configs away from ref state, example is to large
#to do by hand
#probably something to do with building p spaces and q spaces


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
                              thresh_foi=1e-7,    # Threshold for keeping terms when defining FOIS
                              thresh_asci=0.001,     # Threshold of P-space configs to search from
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
    @printf("TCI Thresh: %8.6f  Dim:%8d\n",thresh_cipsi,size(v0)[1])
    #@printf("TCI Thresh: %8.6f  Dim:%8d\n",thresh_cipsi,size(v0)[1])
    println()
    @printf("TCI %5s %12s %12s\n", "Root", "E(0)", "E(2)") 
    for r in 1:nroots
        @printf("TCI %5s %12.8f %12.8f\n",r, e0[r] + ecore, e0[r] + e2[r] + ecore)
    end
end





