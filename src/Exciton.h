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
		// Derived object classes must have a static constant string variable that will be set to a 
		// unique name for identifying different object types when given only an Object base class pointer.
		static const std::string object_type;

		// Derived object constructors must call the base class constructor as well.
		Exciton(const double time, const int tag_num, const KMC_Lattice::Coords& start_coords) : Object(time, tag_num, start_coords) {}
		
		// Derived object classes must define the getObjectType function to retrive the static object type string.
		std::string getObjectType() const { return object_type; }

		// -----------------------------------------------------------------------------------------------
		// Object event classes - One should declare all derived event classes for each type of event
		// that the derived object can perform within the derived object class with public scope.
		// -----------------------------------------------------------------------------------------------

		// This derived event class represents an an exciton creation event
		class Creation : public KMC_Lattice::Event {
		public:
			// Derived event classes must have a static constant string variable that will be set to a 
			// unique name for identifying different event types when given only an Event base class pointer.
			static const std::string event_type;

			// Constructs an empty event that is uninitialized. 
			// Derived event constructors must call the base class constructor.
			Creation() : KMC_Lattice::Event() {}

			// Constructs and initializes an exciton creation event. 
			// Derived event constructors must call the base class constructor.
			Creation(KMC_Lattice::Simulation* sim_ptr) : KMC_Lattice::Event(sim_ptr) {}

			// Derived event classes must define the getEventType function to retrive the static event type string.
			std::string getEventType() const { return event_type; }
		private:

		};

		// This derived event class represnts an exciton hop (move) event
		class Hop : public KMC_Lattice::Event {
		public:
			// Derived event classes must have a static constant string variable that will be set to a 
			// unique name for identifying different event types when given only an Event base class pointer.
			static const std::string event_type;

			// Constructs and initializes an exciton hop event. 
			// Derived event constructors must call the base class constructor.
			Hop(KMC_Lattice::Simulation* sim_ptr) : KMC_Lattice::Event(sim_ptr) {}

			// Derived event classes can have custom calculateRateConstant functions that allows users
			// to implement complex rate equations for the event
			// Here this exciton hop event implements a Forster resonant energy transfer mechanism rate equation.
			void calculateRateConstant(const double prefactor, const double distance, const double E_delta) {
				rate_constant = prefactor * KMC_Lattice::intpow(1.0 / distance, 6);
				if (E_delta > 0) {
					rate_constant *= exp(-E_delta / (KMC_Lattice::K_b*sim_ptr->getTemp()));
				}
			}

			// Derived event classes must define the getEventType function to retrive the static event type string.
			std::string getEventType() const { return event_type; }
		private:

		};

		// This derived event class represents an exciton recombination (disappearance) event
		class Recombination : public KMC_Lattice::Event {
		public:
			// Derived event classes must have a static constant string variable that will be set to a 
			// unique name for identifying different event types when given only an Event base class pointer.
			static const std::string event_type;

			// Constructs and initializes an exciton hop event. 
			// Derived event constructors must call the base class constructor.
			Recombination(KMC_Lattice::Simulation* simulation_ptr) : Event(simulation_ptr) {}

			// Derived event classes must define the getEventType function to retrive the static event type string.
			std::string getEventType() const { return event_type; }

		private:
		};

	private:
	};

}

#endif // EXCITON_H
