# PersistentHomologyComputation
This repo contains source codes for the computation of persistent homology

External dependencies:
- OpenCV 3.3: compile from source (https://opencv.org/releases/). New versions are not tested.
- blitz: download link (https://github.com/blitzpp/blitz). For your convenience, we provide the version we use in this repo.
- Anaconda3: needed for calling c++ routines from python.

To use this library from c++, we use SWIG (https://www.swig.org/) as the wrapper. The interface of this library is defined in "PersistenceComputer.h". The python wrapper is "PersistenceComputer.i". This library is compiled and tested with a Windows 10 system.
