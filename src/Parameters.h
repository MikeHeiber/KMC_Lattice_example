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

	// This derived parameters class extends the Parameters_Simulation base class to store
	// all parameters needed by the simulation.
	class Parameters : public KMC_Lattice::Parameters_Simulation {

	public:
		// Default constructor creates a Parameters object with default values
		Parameters();

		// This function checks whether the current parameter values are valid. Users can expand this
		// function to check for a wide variety of parameter conflicts and to restrict
		// the range of values that are acceptable for users to use as simulation inputs.
		bool checkParameters() const;

		// This function parses the parameters from a parameter file.  Users can design this function
		// to parse whatever type of parameter file format they prefer to use.
		bool importParameters(std::ifstream& inputfile);

		// -----------------------------------------------------------------------------------------------
		// Test Parameters - Users should define what types of tests can be run using the simulation and
		// create boolean "enable" parameters that users will use to select one of the possible tests.
		// -----------------------------------------------------------------------------------------------

		// Here, we only implement one test, an exciton diffusion test.
		// This boolean parameter enables or disables the exciton diffusion test
		bool Enable_diffusion_test = true;

		// The N_tests is a flexible counter for controlling the length of a given simulation test.
		// Here, in the exciton diffusion test we use it to define how many excitons will be tested.
		int N_tests = 0;

		// -----------------------------------------------------------------------------------------------
		// Object Parameters - Users should define the parameters that represent the properties of each
		// type of object in the simulation.
		// -----------------------------------------------------------------------------------------------
		
		// This parameter defines the rate of exciton generation
		double Exciton_generation_rate = 0.0; // (cm^-3 s^-1)
		// This parameter defines the average lifetime of the excitons
		double Exciton_lifetime = 0.0; // (s)
		// This parameter defines the exciton hop rate prefactor (attempt-to-hop frequency)
		double R_exciton_hopping = 0.0; // (s^-1)
		// This parameter defines the exciton hop range cutoff
		int FRET_cutoff = 0; // (nm)

		// -----------------------------------------------------------------------------------------------
		// Lattice Parameters - Users should define the parameters that represent nay added properties of 
		// the lattice sites
		// -----------------------------------------------------------------------------------------------

		// This parameter is used to enable an uncorrelated Gaussian density of states model for the site
		// energies
		bool Enable_gaussian_dos = false;
		// This parameter defines the standard deviation of the Gaussian distribution
		double Site_energy_stdev = 0.0; // (eV)
		// This parameters is used to enable a exponential tail density of states model for the site 
		// energies
		bool Enable_exponential_dos = false;
		// This parameter defines the shape of exponential tail using the so-called Urbach energy.
		double Site_energy_urbach = 0.0; // (eV)

	private:

	};

}

#endif // PARAMETERS_H
