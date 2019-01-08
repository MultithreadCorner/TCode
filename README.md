# TCode

## What is it?
TCode is a C++14 compliant application to simulate the response of solid state sensors in massively parallel platforms on Linux systems. TCode is implemented on top of [Hydra](https://github.com/MultithreadCorner/Hydra) and as such, it can run on OpenMP, CUDA and TBB compatible devices. 
TCode is still in its alpha version and the repository it is taking shape, for the moment we focused on raw performance and will make more complete and user friendly in the near future.

## How it works
TCode uses external 3D maps of electric fields, carrier mobilities and weighting field and energy deposit to simulate the response in current of solid state sensors. The motion of the individual carriers produced in the initial deposit is determined using a 4th other Runge-Kutta using the electric field and the mobilities and assuming that the carriers always move at drift velocity. At each time interval the current induced in the electrod is calculated from the carriers velocity using the corresponding weighting field, according to the [Shockley](https://aip.scitation.org/doi/10.1063/1.1710367)-[Ramo](https://ieeexplore.ieee.org/document/1686997) theorem.

## Dependencies
TCode depends on [ROOT >= v.6.14](https://github.com/root-project/root), [libconfig >= v1.5](https://hyperrealm.github.io/libconfig/) and [TCLAP >= v1.2.1](http://tclap.sourceforge.net/).

## Manual
IN PREPARATION

## Authors
TCode was created by [Andrea Contu](https://github.com/acontu) and [Angelo Loi](https://github.com/angeloloi19).

## Acknowledgement
TCode was developed within the Timespot collaboration, supported by the [Instituto Nazionale di Fisica Nucleare (INFN)](http://home.infn.it/en/).
