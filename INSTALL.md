Build application
-------------------------------

1. clone the git repository: `git clone https://github.com/MultithreadCorner/TCode.git`
2. go to TCode directory: `cd TCode`
3. create a build directory: `mkdir build` 
4. cd build
5. `cmake -DHYDRA_INCLUDE_DIR=<path-to-hydra>../`
6. `make -j8`

Separate executables for CPP, TBB, OMP and CUDA backend will be produced (if available in your machine).
