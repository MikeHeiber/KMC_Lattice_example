// Copyright (c) 2017-2019 Michael C. Heiber
// This source file is part of the KMC_Lattice_example project, which is subject to the MIT License.
// For more information, see the LICENSE file that accompanies this software.
// The KMC_Lattice_example project can be found on Github at https://github.com/MikeHeiber/KMC_Lattice_example

#include "Exciton.h"

using namespace std;

namespace KMC_Lattice_example {

	// Static object_type and event_type strings must be initialized in the .cpp file.
	const string Exciton::object_type = "Exciton";
	const string Exciton::Creation::event_type = "Exciton Creation";
	const string Exciton::Hop::event_type = "Exciton Hop";
	const string Exciton::Recombination::event_type = "Exciton Recombination";

}
