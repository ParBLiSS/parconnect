# parConnect #

### Authors ###

* Chirag Jain

### Install ###


The repository and external submodules can be cloned directly:

    git clone --recursive <GITHUB_URL>
    mkdir build_directory && cd build_directory
    cmake ../parconnect
    make -j4

### Run ###

Inside the build directory, 

    mpirun -np <COUNT OF PROCESSES> ./bin/<EXECUTABLE> <ARGUMENTS>
