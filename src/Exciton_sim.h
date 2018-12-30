// Copyright (c) 2017-2019 Michael C. Heiber
// This source file is part of the KMC_Lattice_example project, which is subject to the MIT License.
// For more information, see the LICENSE file that accompanies this software.
// The KMC_Lattice_example project can be found on Github at https://github.com/MikeHeiber/KMC_Lattice_example

#ifndef EXCITON_SIM_H
#define EXCITON_SIM_H

#include "Event.h"
#include "Exciton.h"
#include "Object.h"
#include "Parameters.h"
#include "Simulation.h"
#include "Utils.h"

namespace KMC_Lattice_example {

	class Site_OSC : public KMC_Lattice::Site {
	public:
		void setEnergy(const double energy_in) { energy = energy_in; }
		double getEnergy() const { return energy; }
	private:
		double energy = 0.0;
	};

	class Exciton_sim : public KMC_Lattice::Simulation {
	public:
		// Functions
		Exciton_sim(const Parameters& params, const int id);
		double calculateDiffusionLength_avg();
		double calculateDiffusionLength_stdev();
		bool checkFinished() const;
		bool executeNextEvent();
		std::vector<double> getDiffusionData();
		int getN_excitons_created();
		int getN_excitons_recombined();
		void outputStatus() const;
	protected:
		std::vector<Site_OSC> sites;
		std::list<Exciton> excitons;
		std::list<Exciton_Hop> exciton_hop_events;
		std::list<Exciton_Recombination> exciton_recombination_events;
	private:
		Parameters params;
		// Derived Paramters
		double R_exciton_generation;
		// Additional Data Structures
		std::vector<double> diffusion_distances;
		Exciton_Creation exciton_creation_event;
		std::list<KMC_Lattice::Event*>::const_iterator exciton_creation_it;
		// Additional Counters
		int N_excitons_created;
		int N_excitons_recombined;
		int N_excitons;
		// Additional Functions
		KMC_Lattice::Coords calculateExcitonCreationCoords();
		void calculateExcitonEvents(Exciton* exciton_it);
		void deleteExciton(Exciton* exciton_ptr);
		bool executeExcitonCreation(const std::list<KMC_Lattice::Event*>::const_iterator event_it);
		bool executeExcitonHop(const std::list<KMC_Lattice::Event*>::const_iterator event_it);
		bool executeExcitonRecombination(const std::list<KMC_Lattice::Event*>::const_iterator event_it);
		std::list<Exciton>::iterator getExcitonIt(const Exciton* exciton_ptr);
		double getSiteEnergy(const KMC_Lattice::Coords& coords) const;
	};

}

#endif // EXCITON_SIM_H
