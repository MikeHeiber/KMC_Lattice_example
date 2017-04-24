#include "Exciton.h"

// Initialize static class members
const string Exciton::name = "Exciton";
const string Exciton_Creation::name = "Exciton Creation";
const string Exciton_Hop::name = "Exciton Hop";
const string Exciton_Recombine::name = "Exciton Recombine";
double Exciton::R_hop = 0;
double Exciton::lifetime = 0;
int Exciton::FRET_cutoff = 0;
