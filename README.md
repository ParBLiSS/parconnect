Parallel Implementation for Computing Graph Connectivity
========================================================================
[![Apache 2.0 License](https://img.shields.io/badge/license-Apache%20v2.0-blue.svg)](LICENSE)

This library implements a distributed connectivity algorithm for large graphs. It can compute the connected components in the undirected graphs or [weakly connected components](http://mathworld.wolfram.com/WeaklyConnectedComponent.html) in the directed graphs. The algorithm is implemented in `C++11` and `MPI`. The codebase supports the graph generation for building de Bruijn graphs from DNA sequence files, synthetic kronecker graphs and also a parallel edgelist file reader for any generic graph.

The coloring algorithm implemented by this codebase is described in the following peer-reviewed publication. Please cite this paper, when using our code for academic purposes:
> **Patrick Flick, Chirag Jain, Tony Pan and Srinivas Aluru.** "A parallel connectivity algorithm for de Bruijn graphs in metagenomic applications." *Proceedings of the International Conference for High Performance Computing, Networking, Storage and Analysis*, ACM, 2015. [![dx.doi.org/10.1145/2807591.2807619](https://img.shields.io/badge/doi-10.1145%2F2807591.2807619-blue.svg)](http://dx.doi.org/10.1145/2807591.2807619)

## Code organization

### Dependencies

- `cmake` version >= 2.6
- `g++` (version 4.8.1+)
- an `MPI` implementation supporting `MPI-2` or `MPI-3`.
- external (third-party) softwares are included in the [`ext/`](`ext/`) directory of this project.

### Download and compile


The repository and external submodules can be downloaded using the recursive clone.

```sh
git clone --recursive <GITHUB_URL>
```

Compile using the cmake tool.

```sh
mkdir build_directory && cd build_directory
cmake ../parconnect
make -j4
```

### Run

Inside the build directory, 

```sh
mpirun -np <COUNT OF PROCESSES> ./bin/<EXECUTABLE> <ARGUMENTS>
```
