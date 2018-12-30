// Copyright (c) 2017-2019 Michael C. Heiber
// This source file is part of the KMC_Lattice_example project, which is subject to the MIT License.
// For more information, see the LICENSE file that accompanies this software.
// The KMC_Lattice_example project can be found on Github at https://github.com/MikeHeiber/KMC_Lattice_example

#ifndef EXCITON_H
#define EXCITON_H

#include "Utils.h"
#include "Object.h"
#include "Event.h"
#include <string>

namespace KMC_Lattice_example {

	class Exciton : public KMC_Lattice::Object {
	public:
		static const std::string object_type;
		Exciton(const double time, const int tag_num, const KMC_Lattice::Coords& start_coords) : Object(time, tag_num, start_coords) {}
		std::string getObjectType() const { return object_type; }
	private:
	};

	class Exciton_Creation : public KMC_Lattice::Event {
	public:
		static const std::string event_type;
		// Constructs an empty event that is uninitialized.
		Exciton_Creation() : KMC_Lattice::Event() {}
		// Constructs and initializes an exciton creation event
		Exciton_Creation(KMC_Lattice::Simulation* sim_ptr) : KMC_Lattice::Event(sim_ptr) {}
		// Gets the event type string that denoes which derived event class it is
		std::string getEventType() const { return event_type; }
	private:

	};

	class Exciton_Hop : public KMC_Lattice::Event {
	public:
		static const std::string event_type;
		Exciton_Hop(KMC_Lattice::Simulation* sim_ptr) : KMC_Lattice::Event(sim_ptr) {}
		void calculateRateConstant(const double prefactor, const double distance, const double E_delta) {
			rate_constant = prefactor * KMC_Lattice::intpow(1.0 / distance, 6);
			if (E_delta > 0) {
				rate_constant *= exp(-E_delta / (KMC_Lattice::K_b*sim_ptr->getTemp()));
			}
		}
		std::string getEventType() const { return event_type; }
	private:

	};

	class Exciton_Recombination : public KMC_Lattice::Event {
	public:
		static const std::string event_type;

		Exciton_Recombination(KMC_Lattice::Simulation* simulation_ptr) : Event(simulation_ptr) {}

		//! \brief Gets the event type string that denotes what type of Event class this is.
		//! \returns The string "Exciton_Recombination".
		std::string getEventType() const { return event_type; }

	private:
	};

}

#endif // EXCITON_H
