// Copyright (c) 2017-2019 Michael C. Heiber
// This source file is part of the KMC_Lattice_example project, which is subject to the MIT License.
// For more information, see the LICENSE file that accompanies this software.
// The KMC_Lattice_example project can be found on Github at https://github.com/MikeHeiber/KMC_Lattice_example

#ifndef PARAMETERS_H
#define PARAMETERS_H

#include "Parameters_Simulation.h"
#include "Utils.h"
#include <fstream>
#include <vector>
#include <string>

namespace KMC_Lattice_example {

	class Parameters : public KMC_Lattice::Parameters_Simulation {

	public:
		// Functions
		Parameters();
		virtual ~Parameters();
		bool checkParameters() const;
		bool importParameters(std::ifstream& inputfile);

		// Variables
		// Test Parameters
		bool Enable_diffusion_test;
		int N_tests;
		// Object Parameters
		// Excitons
		double Exciton_generation_rate; // 1/seconds
		double Exciton_lifetime; // seconds
		double R_exciton_hopping;
		int FRET_cutoff; // nm
		// Lattice Site Parameters
		bool Enable_gaussian_dos;
		double Site_energy_stdev; // eV
		bool Enable_exponential_dos;
		double Site_energy_urbach;

	private:

	};

}

#endif // PARAMETERS_H
