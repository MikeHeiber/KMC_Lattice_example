// Copyright (c) 2017 Michael C. Heiber
// This source file is part of the KMC_Lattice_example project, which is subject to the MIT License.
// For more information, see the LICENSE file that accompanies this software.
// The KMC_Lattice_example project can be found on Github at https://github.com/MikeHeiber/KMC_Lattice_example

#include "Exciton.h"

// Initialize static class members
const string Exciton::name = "Exciton";
const string Exciton_Creation::name = "Exciton Creation";
const string Exciton_Hop::name = "Exciton Hop";
const string Exciton_Recombine::name = "Exciton Recombine";
double Exciton::R_hop = 0;
double Exciton::lifetime = 0;
int Exciton::FRET_cutoff = 0;
