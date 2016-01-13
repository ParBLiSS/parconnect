# parConnect #

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
