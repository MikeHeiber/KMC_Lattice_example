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

	// Derived Site class that adds a site energy property to the lattice.
	// For simple additions to the Site base class, as is the case here, the derived class can 
	// be completely defined quickly in the simulation class header. 
	class Site_OSC : public KMC_Lattice::Site {
	public:
		void setEnergy(const double energy_in) { energy = energy_in; }
		double getEnergy() const { return energy; }
	private:
		double energy = 0.0;
	};

	// Derived Simulation class for simulating exciton creation, diffusio, and recombination
	class Exciton_sim : public KMC_Lattice::Simulation {
	public:
		// Constructor that creates an Exciton_sim object with the specified input parameters
		Exciton_sim(const Parameters& params, const int id);

		// Calculates the average exciton diffusion length of all excitons that have been created 
		// and recombined so far
		double calculateDiffusionLength_avg();

		// Calculates the standard deviation of the exciton diffusion length of all excitons that 
		// have been created and recombined so far
		double calculateDiffusionLength_stdev();

		// This function must be defined in any derived Simulation class to determine when the 
		// simulation test is complete
		bool checkFinished() const;

		// This function is designed to be called from main to execute one iteration of the KMC 
		// algorithm
		bool executeNextEvent();

		// Get a copy of the diffusion length data vector with the displacement distance of all 
		// excitons that have been created and recombined so far
		std::vector<double> getDiffusionData();

		// Gets the number of excitons that have been created so far
		int getN_excitons_created();

		// Gets the number of excitons that have recombined so far
		int getN_excitons_recombined();

		// Outputs the current status of the simulation to the command line
		void outputStatus() const;

	protected:
		// -----------------------------------------------------------------------------------------------
		// Site storage - One needs to store all sites that make up the lattice in the 
		// derived simulation class. Pointers to these sites will be used by the Simulation
		// base class. Use vector because content is constant during the simulation.
		// -----------------------------------------------------------------------------------------------

		// Vector of sites made up of the derived Site_OSC class objects
		std::vector<Site_OSC> sites;

		// -----------------------------------------------------------------------------------------------
		// Object storage - One needs to store each type of object in the simulation. 
		// Pointers to these objects will by used by the Simulation base class.
		// Use list for more rapid adding and removing of objects during the simulation. 
		// -----------------------------------------------------------------------------------------------

		// List of all Exciton objects currently in simulation. 
		std::list<Exciton> excitons;

		// -----------------------------------------------------------------------------------------------
		// Event storage - One needs to store each type of event in the simulation.
		// Pointers to these events will be used bt the Simulation base class.
		// Use list for more rapid adding and removing of events during the simulation.
		// -----------------------------------------------------------------------------------------------

		// List of all Exciton_Hop events in the simulation, one for each Exciton currently in the 
		// simulation.
		std::list<Exciton::Hop> exciton_hop_events;

		// List of all Exiton_Recombunation events in the simulation, one for each Exciton currently 
		// in the simulation.
		std::list<Exciton::Recombination> exciton_recombination_events;

		// Single Exciton_Creation event - No list is needed because there is only ever one of these 
		// events in the simulation.
		Exciton::Creation exciton_creation_event;

	private:
		// Use the derived Parameters_Simulation class object to store all of the input parameters
		Parameters params;

		// -----------------------------------------------------------------------------------------------
		// Derived Parameters - One can also define derived parameters to make things easier.
		// Derived parameters should normally be initialized within the simulation class constructor
		// -----------------------------------------------------------------------------------------------

		// Defines the rate contant for exciton generation, which is calculated based on the input params
		double R_exciton_generation;

		// -----------------------------------------------------------------------------------------------
		// Additional Data Structures - One can define a variety of additional data structures for storing 
		// data needed by any of the simulation tests.
		// -----------------------------------------------------------------------------------------------

		// In this simulation, there is only one test, called the exiton diffusion test
		// To generate diffusion length results, here we use a simple vector for storing the displacement
		// distance of each Exciton once it recombines
		std::vector<double> diffusion_distances;

		// -----------------------------------------------------------------------------------------------
		// Additional Counters - One can define a variety of additional counters to keep track of how many 
		// of each object is on the simulation and how often various events occur during the simulation.
		// -----------------------------------------------------------------------------------------------

		// Keep track of how many Exciton objects are currently in the simulation.
		int N_excitons = 0;

		// Keep track of how many Exciton objects have been created in the simulation so far.
		int N_excitons_created = 0;

		// Keep track of how many Exciton_Recombination events have occurd so far.
		int N_excitons_recombined = 0;

		// -----------------------------------------------------------------------------------------------
		// Calculate event functions - One should define "calculate events" functions for each type of
		// object in the simulation that will calculate all of the possible events for each object type
		// -----------------------------------------------------------------------------------------------

		// Calculates all possible events for the specified Exciton that are declared in the Exciton class. 
		void calculateExcitonEvents(Exciton* exciton_it);

		// -----------------------------------------------------------------------------------------------
		// Execute event functions - One should define "execute event" functions for each type of event
		// for each type of Object in the simulaiton.
		// -----------------------------------------------------------------------------------------------

		// Execute the exciton creation event and return true if the event is successful or false an 
		// error occurs.
		bool executeExcitonCreation(const std::list<KMC_Lattice::Event*>::const_iterator event_it);

		// Execute the exciton hop event and return true if the event is successful or false an error 
		// occurs.
		bool executeExcitonHop(const std::list<KMC_Lattice::Event*>::const_iterator event_it);

		// Execute the exciton recombination event and return true if the event is successful or false an 
		// error occurs.
		bool executeExcitonRecombination(const std::list<KMC_Lattice::Event*>::const_iterator event_it);

		// -----------------------------------------------------------------------------------------------
		// Delete object functions - One should define "delete object" functions for each object type in 
		// the simulation so that "execute event" functions can call them when a particular object needs 
		// to be removed from the simulation.
		// -----------------------------------------------------------------------------------------------

		// This function deletes an Exciton from the simulation
		// The code in this function could be rolled into the executeExcitonRecombination if desired here.
		// But in more complex simulations multiple event may cause excitons to be removed from the simulation,
		// and so resuable delete functions are good to have.
		void deleteExciton(Exciton* exciton_ptr);

		// -----------------------------------------------------------------------------------------------
		// Additional Functions - One should define several types of private functions to make program
		// flow easier to understand and manage.
		// -----------------------------------------------------------------------------------------------

		// This function randomly determines the coordinates where a new exciton will be created
		// The code in this function could be rolled into the executeExcitonCreation function if desired.
		KMC_Lattice::Coords calculateExcitonCreationCoords();

		// This utility function provide an easy resuseable conversion from Exiton pointer to Exciton 
		// list iterator for an alteranteve way to access the Exciton object
		std::list<Exciton>::iterator getExcitonIt(const Exciton* exciton_ptr);

		// This utiliy function provides an easier reuseable way to get the energy of the site at the
		// specified coordinates.
		double getSiteEnergy(const KMC_Lattice::Coords& coords) const;
	};

}

#endif // EXCITON_SIM_H
