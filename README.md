# KMC_Lattice_example v1.0
This object-oriented C++ software tool provides a simple demonstration of how to use the [KMC_Lattice package](https://github.com/MikeHeiber/KMC_Lattice) to create a lattice kinetic Monte Carlo simulation.  This example shows how to extend the base classes in the KMC_Lattice package to simulate exciton creation, diffusion, and decay in an organic semiconducting material.  One could also use this as a simple random-walk simulation.

### Compiling
Compiling requires the boost library for random number generation and an MPI library for parallel processing.

More information about these packages can be found here:
- http://www.boost.org/
- http://www.mpich.org/, http://www.open-mpi.org/, http://mvapich.cse.ohio-state.edu/

### Usage
KMC_Lattice_example.exe takes one required input argument, which is the filename of the input parameter file.

An example parameter file is provided with parameters_default.txt

As an example, to create a single simulation instance on a single processor, the command is:
>    KMC_Lattice_example.exe parameters_default.txt

To run in a parallel processing environment and create 10 simulations on 10 processors to gather more statistics, an example run command is:
>    mpiexec -n 10 KMC_Lattice_example.exe parameters_default.txt

MPI execution commands can be implemented into batch scripts for running KMC_Lattice_example in a supercomputing environment.

### Output
KMC_Lattice_example will create several output files:
- results#.txt -- This text file will contain the results for each processor where the # will be replaced by the processor ID.
- analysis_summary.txt -- When using MPI, this text file will contain average final results from all of the processors.
