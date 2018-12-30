// Copyright (c) 2017-2019 Michael C. Heiber
// This source file is part of the KMC_Lattice_example project, which is subject to the MIT License.
// For more information, see the LICENSE file that accompanies this software.
// The KMC_Lattice_example project can be found on Github at https://github.com/MikeHeiber/KMC_Lattice_example

#include "Exciton_sim.h"

using namespace std;
using namespace KMC_Lattice;

namespace KMC_Lattice_example {

	Exciton_sim::Exciton_sim(const Parameters& params_in, const int id) {
		// Check validity of input parameters
		if (!params_in.checkParameters()) {
			throw invalid_argument("Error! Cannot create Exciton_sim object because the input parameters are invalid.");
		}
		// Set parameters
		params = params_in;
		// Set parameters of Simulation base class using the init function
		init(params, id);
		// Initialize counters
		N_excitons = 0;
		N_excitons_created = 0;
		N_excitons_recombined = 0;
		// Initialize lattice sites
		Site_OSC site;
		site.clearOccupancy();
		sites.assign(lattice.getNumSites(), site);
		// Initialize site energies using Utils functions
		vector<double> site_energies(lattice.getNumSites());
		if (params.Enable_gaussian_dos) {
			createGaussianDOSVector(site_energies, 0, params.Site_energy_stdev, generator);
		}
		else if (params.Enable_exponential_dos) {
			createExponentialDOSVector(site_energies, 0, params.Site_energy_urbach, generator);
		}
		// Set the site energies
		for (int i = 0; i < lattice.getNumSites(); i++) {
			if (params.Enable_gaussian_dos || params.Enable_exponential_dos) {
				sites[i].setEnergy(site_energies[i]);
			}
		}
		// Set the site pointers for the lattice
		vector<Site*> site_ptrs((int)sites.size());
		for (int i = 0; i < (int)sites.size(); i++) {
			site_ptrs[i] = &sites[i];
		}
		lattice.setSitePointers(site_ptrs);
		// Initialize exciton creation event
		R_exciton_generation = params.Exciton_generation_rate * lattice.getNumSites()*intpow(1e-7*lattice.getUnitSize(), 3);
		exciton_creation_event = Exciton_Creation(this);
		exciton_creation_event.calculateExecutionTime(R_exciton_generation);
		exciton_creation_it = addEvent(&exciton_creation_event);
	}

	double Exciton_sim::calculateDiffusionLength_avg() {
		return vector_avg(diffusion_distances);
	}

	double Exciton_sim::calculateDiffusionLength_stdev() {
		return vector_stdev(diffusion_distances);
	}

	Coords Exciton_sim::calculateExcitonCreationCoords() {
		// Faster method of choosing a random site.  This becomes slow when the lattice has high occupancy.
		// Use faster method when lattice is less than 50% occupied.
		int N_tries = 0;
		while (N_excitons < (0.5*lattice.getLength()*lattice.getWidth()*lattice.getHeight()) && N_tries < 10) {
			auto dest_coords = lattice.generateRandomCoords();
			if (!lattice.isOccupied(dest_coords)) {
				return dest_coords;
			}
			N_tries++;
		}
		// Method of choosing one of the empty sites.  This is slowish becuase it must loop through all sites first
		// Get vector of possible creation sites
		vector<long int> indices;
		indices.reserve(lattice.getLength()*lattice.getWidth()*lattice.getHeight());
		for (int x = 0; x < lattice.getLength(); x++) {
			for (int y = 0; y < lattice.getWidth(); y++) {
				for (int z = 0; z < lattice.getHeight(); z++) {
					Coords dest_coords(x, y, z);
					if (!lattice.isOccupied(dest_coords)) {
						indices.push_back(lattice.getSiteIndex(dest_coords));
					}
				}
			}
		}
		// Print an error message if no unoccupied sites are found
		if ((int)indices.size() == 0) {
			cout << getId() << ": Error! An empty site for exciton creation could not be found." << endl;
			return Coords(-1, -1, -1);
		}
		uniform_int_distribution<long int> distn(0, (long int)indices.size());
		return lattice.getSiteCoords(indices[distn(generator)]);
	}

	// Each object type should have an event calculation function to that makes sure all possible event transitions are calculated
	void Exciton_sim::calculateExcitonEvents(Exciton* exciton_ptr) {
		// Gather information about the exciton
		const auto exciton_it = getExcitonIt(exciton_ptr);
		const Coords object_coords = exciton_it->getCoords();
		// Create vector to store pointers to all possible events that the input exciton can do
		vector<Event*> possible_events;
		// Calculate all possible Exciton Hopping events and add their pointers to the possible_events vector
		// Assess the nearby sites to determine if a hop can occur to them and calculates what rate constant is to each
		static Exciton_Hop event_hop(this);
		// Exciton Hop range is calculated in lattice units based on the specified hop cutoff distance in real space units
		static const int range = (int)ceil((double)params.FRET_cutoff / lattice.getUnitSize());
		static const int dim = (2 * range + 1);
		// Use pre-allocated event vector so that all events in the search do not need to be created each time
		static vector<Exciton_Hop> hops_temp(dim*dim*dim, event_hop);
		for (int i = -range; i <= range; i++) {
			for (int j = -range; j <= range; j++) {
				for (int k = -range; k <= range; k++) {
					// Use the Lattice class checkMoveValidity function to see if a move with displacement (i,j,k) is possible
					if (!lattice.checkMoveValidity(object_coords, i, j, k)) {
						continue;
					}
					// Use the Lattice class calculateDestinationCoords functions to determine the destination coordinates of the proposed move
					Coords dest_coords;
					lattice.calculateDestinationCoords(object_coords, i, j, k, dest_coords);
					// Use the Lattice class isOccupied function to check if the site at the destination coordinates is unoccupied
					if (lattice.isOccupied(dest_coords)) {
						continue;
					}
					// Calculate the real space distance of the proposed move
					double distance = lattice.getUnitSize()*sqrt((double)(i*i + j * j + k * k));
					// Save the move as a possible event only if the move distance is less than the specified cutoff distance
					if (!((distance - 0.0001) > params.FRET_cutoff)) {
						int index = (i + range)*dim*dim + (j + range)*dim + (k + range);
						// Must specify which object the event is associated with
						hops_temp[index].setObjectPtr(exciton_ptr);
						// Must specify the event destination coords
						hops_temp[index].setDestCoords(dest_coords);
						// Must calculate the event rate constant
						double E_delta = (getSiteEnergy(dest_coords) - getSiteEnergy(object_coords));
						hops_temp[index].calculateRateConstant(params.R_exciton_hopping, distance, E_delta);
						// Save the calculated exciton hop event as a possible event
						possible_events.push_back(&hops_temp[index]);
					}
				}
			}
		}
		// Account for possible Exciton Recombination event
		auto recombination_event_it = exciton_recombination_events.begin();
		advance(recombination_event_it, std::distance(excitons.begin(), exciton_it));
		possible_events.push_back(&*recombination_event_it);
		// Use Simulation class determinePathway function to select which event will be next
		Event* event_ptr_target = determinePathway(possible_events);
		// Copy the chosen temp event to the appropriate main event list and set the target event pointer to the corresponding event from the main list
		string event_type = event_ptr_target->getEventType();
		if (event_type.compare(Exciton_Hop::event_type) == 0) {
			auto hop_list_it = exciton_hop_events.begin();
			std::advance(hop_list_it, std::distance(excitons.begin(), exciton_it));
			*hop_list_it = *static_cast<Exciton_Hop*>(event_ptr_target);
			event_ptr_target = &(*hop_list_it);
		}
		// Set the chosen event for the exciton using the Simulation class setObjectEvent function
		setObjectEvent(exciton_ptr, event_ptr_target);
	}

	// The function must be defined in the derived Simulation class
	bool Exciton_sim::checkFinished() const {
		if (params.Enable_diffusion_test) {
			return (N_excitons_recombined == params.N_tests);
		}
		cout << "Error checking simulation finish conditions." << endl;
		return true;
	}

	void Exciton_sim::deleteExciton(Exciton* exciton_ptr) {
		// Gather exciton information
		auto exciton_it = getExcitonIt(exciton_ptr);
		// remove the Object and Event pointers using Simulation class removeObject function
		removeObject(exciton_ptr);
		// Locate corresponding recombination event
		auto recombination_list_it = exciton_recombination_events.begin();
		advance(recombination_list_it, std::distance(excitons.begin(), exciton_it));
		// Locate corresponding hop event
		auto hop_list_it = exciton_hop_events.begin();
		advance(hop_list_it, std::distance(excitons.begin(), exciton_it));
		// Delete exciton
		excitons.erase(exciton_it);
		// Delete exciton recombination event
		exciton_recombination_events.erase(recombination_list_it);
		// Delete exciton hop event
		exciton_hop_events.erase(hop_list_it);
	}

	// Each event type should have an associated execute function
	bool Exciton_sim::executeExcitonCreation(const list<Event*>::const_iterator event_it) {
		// Determine coords for the new exciton
		const Coords coords_new = calculateExcitonCreationCoords();
		// Create the new exciton and at it to the simulation
		Exciton exciton_new(getTime(), N_excitons_created + 1, coords_new);
		excitons.push_back(exciton_new);
		// Add new exciton to the Simulation class using its addObject function
		addObject(&excitons.back());
		// Add an empty hop and recombine event to the corresponding lists
		Exciton_Hop hop_event(this);
		exciton_hop_events.push_back(hop_event);
		Exciton_Recombination recombination_event(this);
		// Set the recombination event associated object using the Event class setObjectPtr function
		recombination_event.setObjectPtr(&excitons.back());
		// Since the rate constant for all recombination events is the same and does not change,
		// it can be set during initialization using the Event class calculateRateConstant function
		recombination_event.calculateRateConstant(1.0 / params.Exciton_lifetime);
		exciton_recombination_events.push_back(recombination_event);
		// Update exciton counters
		N_excitons_created++;
		N_excitons++;
		// Find all nearby excitons using the Simulation class findRecalcObjects function and calculate their next events
		auto neighbors = findRecalcObjects(coords_new, coords_new);
		for (auto& item : neighbors) {
			calculateExcitonEvents(static_cast<Exciton*>(item));
		}
		// Calculate when the next exciton creation event will occur
		exciton_creation_event.calculateExecutionTime(R_exciton_generation);
		return true;
	}

	// Each event type should have an associated execute function
	bool Exciton_sim::executeExcitonHop(const list<Event*>::const_iterator event_it) {
		if (lattice.isOccupied((*event_it)->getDestCoords())) {
			cout << "Exciton hop cannot be executed. Destination site is already occupied." << endl;
			return false;
		}
		else {
			// Get event info
			Coords coords_initial = ((*event_it)->getObjectPtr())->getCoords();
			// Move the exciton using the Simulation class moveObject function
			moveObject((*event_it)->getObjectPtr(), (*event_it)->getDestCoords());
			// Find all nearby excitons using the Simulation class findRecalcObjects function and calculate their next events
			auto neighbors = findRecalcObjects(coords_initial, (*event_it)->getDestCoords());
			for (auto& item : neighbors) {
				calculateExcitonEvents(static_cast<Exciton*>(item));
			}
			return true;
		}
	}

	// Each event type should have an associated execute function
	bool Exciton_sim::executeExcitonRecombination(const list<Event*>::const_iterator event_it) {
		// Get event info
		int exciton_tag = ((*event_it)->getObjectPtr())->getTag();
		Coords coords_initial = ((*event_it)->getObjectPtr())->getCoords();
		// Output final diffusion displacement distance 
		if (params.Enable_diffusion_test) {
			diffusion_distances.push_back(lattice.getUnitSize()*((*event_it)->getObjectPtr())->calculateDisplacement());
		}
		// Delete Exciton and its events
		deleteExciton(static_cast<Exciton*>((*event_it)->getObjectPtr()));
		// Update exciton counters
		N_excitons--;
		N_excitons_recombined++;
		// Find all nearby excitons using the Simulation class findRecalcObjects function and calculate their next events
		auto neighbors = findRecalcObjects(coords_initial, coords_initial);
		for (auto& item : neighbors) {
			calculateExcitonEvents(static_cast<Exciton*>(item));
		}
		return true;
	}

	bool Exciton_sim::executeNextEvent() {
		// Use the Simulation class chooseNextEvent function to determine which event will be executed
		auto event_it = chooseNextEvent();
		// Gather event info
		string event_name = (*event_it)->getEventType();
		// Update simulation time
		setTime((*event_it)->getExecutionTime());
		// Determine which event type has been chosen and run the appropriate execute function
		if (event_name.compare(Exciton_Creation::event_type) == 0) {
			return executeExcitonCreation(event_it);
		}
		else if (event_name.compare(Exciton_Hop::event_type) == 0) {
			return executeExcitonHop(event_it);
		}
		else if (event_name.compare(Exciton_Recombination::event_type) == 0) {
			return executeExcitonRecombination(event_it);
		}
		else {
			//error
			cout << "Error! Valid event not found when calling executeNextEvent." << endl;
			return false;
		}
	}

	// Simple utility function for converting from object pointer to object list iterator
	list<Exciton>::iterator Exciton_sim::getExcitonIt(const Exciton* exciton_ptr) {
		auto it = find_if(excitons.begin(), excitons.end(), [exciton_ptr](Exciton& a) {return (a.getTag() == exciton_ptr->getTag()); });
		if (it == excitons.end()) {
			cout << "Error! Exciton iterator could not be located." << endl;
		}
		return it;
	}

	vector<double> Exciton_sim::getDiffusionData() {
		return diffusion_distances;
	}

	int Exciton_sim::getN_excitons_created() {
		return N_excitons_created;
	}

	int Exciton_sim::getN_excitons_recombined() {
		return N_excitons_recombined;
	}

	void Exciton_sim::outputStatus() const {
		cout << getId() << ": Time = " << getTime() << " seconds.\n";
		cout << getId() << ": " << N_excitons_created << " excitons have been created and " << getN_events_executed() << " events have been executed.\n";
		cout << getId() << ": There are " << N_excitons << " excitons in the lattice.\n";
		for (auto& item : excitons) {
			cout << getId() << ": Exciton " << item.getTag() << " is at " << item.getCoords().x << "," << item.getCoords().y << "," << item.getCoords().z << ".\n";
		}
		cout.flush();
	}

	// Make it easier to get the energy of a particular site in the lattice
	double Exciton_sim::getSiteEnergy(const Coords& coords) const {
		return sites[lattice.getSiteIndex(coords)].getEnergy();
	}

}
