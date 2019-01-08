Build application
-------------------------------
Make sure that you have the following packages available in your system:
- [Hydra >= v2.2.1](https://github.com/MultithreadCorner/Hydra/releases/tag/v2.2.1)
- [libconfig >= v1.5](https://github.com/hyperrealm/libconfig)
- [TCLAP >= v1.2.1](http://tclap.sourceforge.net/)
- [ROOT >= v6.14](https://root.cern.ch/)

To generate executables to run on nVidia GPUs you also need an installation of [CUDA >= 8.1](https://developer.nvidia.com/cuda-toolkit) and to use a [compatible GCC version](https://docs.nvidia.com/cuda/).

Build TCode following the instructions below:
1. clone the git repository: `git clone https://github.com/MultithreadCorner/TCode.git`
2. go to TCode directory: `cd TCode`
3. create a build directory: `mkdir build` 
4. go to build directory cd `build`
5. run cmake specifying the path to Hydra: `cmake -DHYDRA_INCLUDE_DIR='path-to-hydra'../`
6. compile (all backends): `make -j8`

Several executables, separate for each available backend (CPP, TBB, OMP and CUDA) will be generated.
