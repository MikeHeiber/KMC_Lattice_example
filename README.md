# KMC_Lattice_example
This object-oriented C++11 software tool provides a simple demonstration of how to use v2.0 of the [KMC_Lattice package](https://github.com/MikeHeiber/KMC_Lattice) to create a cubic lattice kinetic Monte Carlo simulation.  This example shows how to extend the base classes in the KMC_Lattice package to simulate exciton creation, diffusion, and decay in an organic semiconducting material.  One could also use this as a simple random-walk simulation.

### Compiling
This software packages makes use of a number of feaures added in C++11, so your compiler must support the C++11 standard.
Compiling also requires an MPI library for parallel processing.

Links to several common MPI packages are provided here:
- http://www.mpich.org/
- http://www.open-mpi.org/
- http://mvapich.cse.ohio-state.edu/

### Usage
KMC_Lattice_example.exe takes one required input argument, which is the filename of the input parameter file.

An example parameter file is provided with parameters_default.txt

To create 10 simulations on 10 processors to gather statistics, an example run command is:
>    mpiexec -n 10 KMC_Lattice_example.exe parameters_default.txt

MPI execution commands can be implemented into batch scripts for running KMC_Lattice_example in a supercomputing environment.

### Output
KMC_Lattice_example will create several output files:
- results#.txt -- This text file will contain the results for each processor where the # will be replaced by the processor ID.
- analysis_summary.txt -- When using MPI, this text file will contain average final results from all of the processors.
