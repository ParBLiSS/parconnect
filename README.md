Parallel Implementation for Computing Graph Connectivity
========================================================
[![Apache 2.0 License](https://img.shields.io/badge/license-Apache%20v2.0-blue.svg)](LICENSE)

### Authors ###

* Chirag Jain

### Install ###


The repository and external submodules can be cloned directly:

```sh
git clone --recursive <GITHUB_URL>
mkdir build_directory && cd build_directory
cmake ../parconnect
make -j4
```

### Run ###

Inside the build directory, 

```sh
mpirun -np <COUNT OF PROCESSES> ./bin/<EXECUTABLE> <ARGUMENTS>
```
