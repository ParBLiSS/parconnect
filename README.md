parConnect
==========

### Install ###


The repository and external submodules can be cloned directly:

    git clone --recursive <GITHUB_URL>
    mkdir build_directory && cd build_directory
    cmake ../parConnect
    make 

### Run ###

Inside the build directory, 

    mpirun -np <COUNT OF PROCESSES> ./bin/<EXECUTABLE> <ARGUMENTS>

### Planned

- Implementation of parallel algorithm for computing connectivity on large graphs.
- Targeted for distributed parallel system.
- Include implementation of Vishkin's PRAM algorithm for connectivity 
- Design algorithm to compute connected components in undirected graphs or weakly connected components in directed graphs. 
- Algorithm should be a mix of BFS and coloring approach.
- Algorithm should perform best on different graph topologies and not just scale-free graphs.
- Compare against boost connectivity algorithm, GraphLab, parallel BFS, coloring. 

