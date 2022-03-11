# Butadiene

1. Run PySCF to get the orbitals

```bash
python rhf.jl
```

This runs both RHF and CIS to generate a state averaged density matrix. 
The partially occupied eigenvectors of this density matrix will be used to define our systems active space. 

2. Generate active space.
```bash
python nat_orb_active_space.py
```
This will generate the orbitals for our active space, and then compute the 1 and 2 electron integrals in this space, 
including the effective 1 particle interactions coming from the frozen doubly occupied orbitals. 
This will also use Pipek-Mezey to localize the active space orbitals. 

3. Cluster active space. 
Next we run `clustering.jl` file to use spectral clustering to automatically group the orbitals into sets
of strongly interacting regions. Here we use the density matrix as an Adjacency matrix, for the graph laplacian. 
There are multiple ways to do this, but ultimately this shouldn't affect the results because it only defines an input guess to the CMF code which optimizes the orbitals. 

```julia
julia --project=./ -tauto 
julia>  include("clustering.jl");
```

using the REPL, or just as a script:

```julia
julia --project=./ clustering.jl >& clustering.out
```
