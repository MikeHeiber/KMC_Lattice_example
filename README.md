<!---
# Copyright (c) 2017-2019 Michael C. Heiber
# This source file is part of the KMC_Lattice_example project, which is subject to the MIT License.
# For more information, see the LICENSE file that accompanies this software.
# The KMC_Lattice_example project can be found on Github at https://github.com/MikeHeiber/KMC_Lattice_example
--->

# KMC_Lattice_example

This object-oriented modern C++ software tool provides a simple demonstration of how to use v2.0 of the [KMC_Lattice library](https://github.com/MikeHeiber/KMC_Lattice) to create a cubic lattice kinetic Monte Carlo (KMC) simulation. 
This example shows how to extend the base classes in the KMC_Lattice library to create derived classes that simulate exciton creation, diffusion, and decay in an organic semiconducting material. 
The source code in this package has detailed comments that explain and show how to implement a KMC simulation using the KMC_Lattice library.
For a more complex and complete example using the KMC_Lattice library to produce computational research software, see the [Excimontec software tool](https://github.com/MikeHeiber/Excimontec).

If you have any questions about the example presented here, please post them in the [Issues](https://github.com/MikeHeiber/KMC_Lattice_example/issues) section. 

#### Building the Executable

This software tool uses [Message Passing Interface (MPI)](https://computing.llnl.gov/tutorials/mpi/) to utilize parallel computing power. 
As a result, using this tool requires that an MPI library is pre-installed on your system, and the final executable must be built on your specific system. 
We cannot provide pre-built binaries for your system. 
Contact your HPC admin to determine the protocols for building MPI applications on your HPC system. 
In many cases, the HPC system will already be configured for you, and the package comes with a default makefile that can be used with the [GCC compiler](https://gcc.gnu.org/) or the [PGI compiler](https://www.pgroup.com/). 

If you wish, you can also install MPI on your own personal workstation and then build the tool there as well. 
More information about common MPI packages can be found here:
- Open MPI, http://www.open-mpi.org/
- MPICH, http://www.mpich.org/
- MVAPICH, http://mvapich.cse.ohio-state.edu/

Once you have an MPI library installed, to build the executable, first copy the project directory to your machine. 
On Linux this can be done using the command,

```git clone --recurse-submodules https://github.com/MikeHeiber/KMC_Lattice_example```

Then set the project directory as your working directory,

```cd KMC_Lattice_example```

and finally build the executable with the default makefile.

```make```

In the default makefile, compilation flags have been set for the GCC and PGI compilers. 
If you are using another compiler, you will need to edit the makefile and define your own compiler options.

Please report any build errors in the [Issues](https://github.com/MikeHeiber/KMC_Lattice_example/issues) section. 

### Usage

KMC_Lattice_example.exe takes one required input argument, which is the filename of the input parameter file.

An example parameter file is provided with parameters_default.txt

To create 10 simulations on 10 processors to gather statistics, an example run command is:
```mpiexec -n 10 KMC_Lattice_example.exe parameters_default.txt```

MPI execution commands can be implemented into batch scripts for running KMC_Lattice_example in a supercomputing environment.

### Output

KMC_Lattice_example will create several output files:
- results#.txt -- This text file will contain the results for each processor where the # will be replaced by the processor ID.
- analysis_summary.txt -- When using MPI, this text file will contain average final results from all of the processors.
