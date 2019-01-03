// Copyright (c) 2017-2019 Michael C. Heiber
// This source file is part of the KMC_Lattice_example project, which is subject to the MIT License.
// For more information, see the LICENSE file that accompanies this software.
// The KMC_Lattice_example project can be found on Github at https://github.com/MikeHeiber/KMC_Lattice_example

#include "Parameters.h"

using namespace std;
using namespace KMC_Lattice;

namespace KMC_Lattice_example {

	Parameters::Parameters() {}

	bool Parameters::checkParameters() const {
		if (!(Temperature > 0)) {
			cout << "Error! The temperature must be greater than zero." << endl;
			return false;
		}
		if (Recalc_cutoff < FRET_cutoff) {
			cout << "Error! The event recalculation cutoff radius must not be less than the FRET cutoff radius." << endl;
			return false;
		}
		if (!(N_tests > 0)) {
			cout << "Error! The number of exciton diffusion tests must be greater than zero." << endl;
			return false;
		}
		if (!(Exciton_generation_rate > 0) || !(Exciton_lifetime > 0) || !(R_exciton_hopping > 0) || !(FRET_cutoff > 0)) {
			cout << "Error! All exciton properties must be greater than zero." << endl;
			return false;
		}
		if (Enable_gaussian_dos && Enable_exponential_dos) {
			cout << "Error! The Gaussian and exponential disorder models cannot both be enabled." << endl;
			return false;
		}
		if (Enable_gaussian_dos && Site_energy_stdev < 0) {
			cout << "Error! When using the Gaussian disorder model, the standard deviation cannot be negative." << endl;
			return false;
		}
		if (Enable_exponential_dos && Site_energy_urbach < 0) {
			cout << "Error! When using the exponential disorder model, the Urbach energy cannot be negative." << endl;
			return false;
		}
		if (Enable_selective_recalc && Recalc_cutoff < FRET_cutoff) {
			cout << "Error! When using the KMC selective recalculation algorithm, the recalculation cutoff distance must not be less than the FRET cutoff distance." << endl;
			return false;
		}
		return true;
	}

	// This parameter import function is designed to parse parameter files based on the format
	// that is present in the parameters_default.txt file included with this package
	bool Parameters::importParameters(ifstream& inputfile) {
		// Read all lines of the file into a string vector
		vector<string> stringvars;
		while (inputfile.good()) {
			string line;
			getline(inputfile, line);
			if ((line.substr(0, 2)).compare("--") != 0 && (line.substr(0, 2)).compare("##") != 0) {
				auto pos = line.find("/", 0);
				auto var = line.substr(0, pos - 1);
				stringvars.push_back(var);
			}
		}
		// Parse each line of the parameter file
		// Boolean "true" or "false" values can be parsed using the Utils str2bool function.
		int i = 0;
		// KMC Algorithm Parameters
		try {
			Enable_FRM = str2bool(stringvars[i]);
		}
		catch (invalid_argument& exception) {
			cout << exception.what() << endl;
			cout << "Error setting first reaction method option." << endl;
			return false;
		}
		i++;
		try {
			Enable_selective_recalc = str2bool(stringvars[i]);
		}
		catch (invalid_argument& exception) {
			cout << exception.what() << endl;
			cout << "Error setting selective recalculation method option." << endl;
			return false;
		}
		i++;
		Recalc_cutoff = atoi(stringvars[i].c_str());
		i++;
		try {
			Enable_full_recalc = str2bool(stringvars[i]);
		}
		catch (invalid_argument& exception) {
			cout << exception.what() << endl;
			cout << "Error setting full recalculation method option." << endl;
			return false;
		}
		i++;
		//enable_periodic_x
		try {
			Params_lattice.Enable_periodic_x = str2bool(stringvars[i]);
		}
		catch (invalid_argument& exception) {
			cout << exception.what() << endl;
			cout << "Error setting x-periodic boundary options" << endl;
			return false;
		}
		i++;
		//enable_periodic_y
		try {
			Params_lattice.Enable_periodic_y = str2bool(stringvars[i]);
		}
		catch (invalid_argument& exception) {
			cout << exception.what() << endl;
			cout << "Error setting y-periodic boundary options" << endl;
			return false;
		}
		i++;
		//enable_periodic_z
		try {
			Params_lattice.Enable_periodic_z = str2bool(stringvars[i]);
		}
		catch (invalid_argument& exception) {
			cout << exception.what() << endl;
			cout << "Error setting z-periodic boundary options" << endl;
			return false;
		}
		i++;
		Params_lattice.Length = atoi(stringvars[i].c_str());
		i++;
		Params_lattice.Width = atoi(stringvars[i].c_str());
		i++;
		Params_lattice.Height = atoi(stringvars[i].c_str());
		i++;
		Params_lattice.Unit_size = atof(stringvars[i].c_str());
		i++;
		Temperature = atoi(stringvars[i].c_str());
		i++;
		Recalc_cutoff = atoi(stringvars[i].c_str());
		i++;
		//Tests
		//enable_exciton_diffusion_test
		try {
			Enable_diffusion_test = str2bool(stringvars[i]);
		}
		catch (invalid_argument& exception) {
			cout << exception.what() << endl;
			cout << "Error setting exciton diffusion test options" << endl;
			return false;
		}
		i++;
		N_tests = atoi(stringvars[i].c_str());
		i++;
		// Exciton Parameters
		Exciton_generation_rate = atof(stringvars[i].c_str());
		i++;
		Exciton_lifetime = atof(stringvars[i].c_str());
		i++;
		R_exciton_hopping = atof(stringvars[i].c_str());
		i++;
		FRET_cutoff = atoi(stringvars[i].c_str());
		i++;
		// Energetic Disorder Parameters
		//enable_gaussian_dos
		try {
			Enable_gaussian_dos = str2bool(stringvars[i]);
		}
		catch (invalid_argument& exception) {
			cout << exception.what() << endl;
			cout << "Error setting Gaussian DOS options" << endl;
			return false;
		}
		i++;
		Site_energy_stdev = atof(stringvars[i].c_str());
		i++;
		//enable_exponential_dos
		try {
			Enable_exponential_dos = str2bool(stringvars[i]);
		}
		catch (invalid_argument& exception) {
			cout << exception.what() << endl;
			cout << "Error setting Exponential DOS options" << endl;
			return false;
		}
		i++;
		Site_energy_urbach = atof(stringvars[i].c_str());
		i++;
		return true;
	}
}
