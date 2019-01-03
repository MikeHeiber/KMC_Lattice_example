# Copyright (c) 2017-2019 Michael C. Heiber
# This source file is part of the KMC_Lattice_example project, which is subject to the MIT License.
# For more information, see the LICENSE file that accompanies this software.
# The KMC_Lattice_example project can be found on Github at https://github.com/MikeHeiber/KMC_Lattice_example

ifeq ($(lastword $(subst /, ,$(CXX))),g++)
	FLAGS += -Wall -Wextra -O3 -std=c++11 -I. -Isrc -IKMC_Lattice/src
endif
ifeq ($(lastword $(subst /, ,$(CXX))),pgc++)
	FLAGS += -O2 -Minform=warn -fastsse -Mvect -std=c++11 -Mdalign -Munroll -Mipa=fast -Kieee -m64 -I. -Isrc -IKMC_Lattice/src
endif

OBJS = src/Exciton_sim.o src/Exciton.o src/Parameters.o

all : KMC_Lattice_example.exe
ifndef FLAGS
	$(error Valid compiler not detected.)
endif

KMC_Lattice_example.exe : src/main.o $(OBJS) KMC_Lattice/libKMC.a
	mpicxx $(FLAGS) $^ -o $@

KMC_Lattice/libKMC.a : KMC_Lattice/src/*.h
	$(MAKE) -C KMC_Lattice

src/main.o : src/main.cpp src/Exciton_sim.h src/Exciton.h src/Parameters.h KMC_Lattice/libKMC.a
	mpicxx $(FLAGS) -c $< -o $@

src/Exciton_sim.o : src/Exciton_sim.cpp src/Exciton_sim.h src/Exciton.h src/Parameters.h KMC_Lattice/libKMC.a
	mpicxx $(FLAGS) -c $< -o $@

src/Parameters.o : src/Parameters.cpp src/Parameters.h KMC_Lattice/libKMC.a
	mpicxx $(FLAGS) -c $< -o $@

src/Exciton.o : src/Exciton.cpp src/Exciton.h KMC_Lattice/libKMC.a
	mpicxx $(FLAGS) -c $< -o $@

clean:
	$(MAKE) -C KMC_Lattice clean
	-rm src/*.o src/*.gcno* src/*.gcda *~ KMC_Lattice_example.exe
