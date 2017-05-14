// Copyright (c) 2017 Michael C. Heiber
// This source file is part of the KMC_Lattice_example project, which is subject to the MIT License.
// For more information, see the LICENSE file that accompanies this software.
// The KMC_Lattice_example project can be found on Github at https://github.com/MikeHeiber/KMC_Lattice_example

#include "Exciton_sim.h"

Exciton_sim::Exciton_sim(const Parameters_Exciton_Sim& params,const int id){
    // Set parameters of Simulation base class
    initializeSimulation(params,id);
    // Set additional parameters
    // Object Parameters
    // Excitons
    Exciton_generation_rate = params.Exciton_generation_rate;
    Exciton_lifetime = params.Exciton_lifetime;
    R_exciton_hopping = params.R_exciton_hopping;
    FRET_cutoff = params.FRET_cutoff;
    // Test Parameters
    Enable_diffusion_test = params.Enable_diffusion_test;
    N_tests = params.N_tests;
    // Lattice Parameters
    Enable_gaussian_dos = params.Enable_gaussian_dos;
    Site_energy_stdev = params.Site_energy_stdev;
    Enable_exponential_dos = params.Enable_exponential_dos;
    Site_energy_urbach = params.Site_energy_urbach;
    // Output files
    //
    // Initialize energetic disorder
    if(Enable_gaussian_dos){
        site_energies.assign(getNumSites(),0);
        createGaussianDOSVector(site_energies,0,Site_energy_stdev,gen);
    }
    else if(Enable_exponential_dos){
        site_energies.assign(getNumSites(),0);
        createExponentialDOSVector(site_energies,0,Site_energy_urbach,gen);
    }
    else{
        site_energies.push_back(0);
    }
    // Initialize counters
    N_excitons = 0;
    N_excitons_created = 0;
    N_excitons_recombined = 0;
    // Initialize sites
    Site_OSC site;
    site.clearOccupancy();
    sites.assign(getNumSites(),site);
    for(int i=0;i<getNumSites();i++){
        if(Enable_gaussian_dos || Enable_exponential_dos){
            sites[i].setEnergyIt(site_energies.begin()+i);
        }
        else{
            sites[i].setEnergyIt(site_energies.begin());
        }
        addSite(&sites[i]);
    }
    // Initialize exciton creation event
    R_exciton_generation = Exciton_generation_rate*getNumSites()*intpow(1e-7*getUnitSize(),3);
    exciton_creation_event.calculateExecutionTime(R_exciton_generation,getTime());
    exciton_creation_it = addEvent(&exciton_creation_event);
}

double Exciton_sim::calculateDiffusionLength_avg(){
    return vector_avg(diffusion_distances);
}

double Exciton_sim::calculateDiffusionLength_stdev(){
    return vector_stdev(diffusion_distances);
}

Coords Exciton_sim::calculateExcitonCreationCoords(){
    Coords dest_coords;
    bool success = false;
    while(!success){
        dest_coords = getRandomCoords();
        if(loggingEnabled()){
            *Logfile << "Attempting to create exciton at " << dest_coords.x << "," << dest_coords.y << "," << dest_coords.z << "." << endl;
        }
        if(isOccupied(dest_coords)){
            continue;
        }
        else{
            success = true;
            return dest_coords;
        }
    }
}

void Exciton_sim::calculateExcitonEvents(const list<Object*>::iterator object_it){
    const auto exciton_it = getExcitonIt(*object_it);
    const Coords object_coords = exciton_it->getCoords();
    if(loggingEnabled()){
        *Logfile << "Calculating events for exciton " << exciton_it->getTag() << " at site " << object_coords.x << "," << object_coords.y << "," << object_coords.z << "." << endl;
    }
    Coords dest_coords;
    double distance,E_delta,rate;
    int index;
    // Exciton Hopping
    static Exciton_Hop event_hop;
    static const int range = ceil(FRET_cutoff/getUnitSize());
    static const int dim = (2*range+1);
    static vector<Exciton_Hop> hops_temp(dim*dim*dim,event_hop);
    static vector<bool> hops_valid(dim*dim*dim,false);
    for(int i=-range;i<=range;i++){
        for(int j=-range;j<=range;j++){
            for(int k=-range;k<=range;k++){
                index = (i+range)*dim*dim+(j+range)*dim+(k+range);
                if(!checkMoveEventValidity(object_coords,i,j,k)){
                    hops_valid[index] = false;
                    continue;
                }
                dest_coords = calculateDestinationCoords(object_coords,i,j,k);
                if((*getSiteIt(dest_coords))->isOccupied()){
                    hops_valid[index] = false;
                    continue;
                }
                distance = getUnitSize()*sqrt((double)(i*i+j*j+k*k));
                if(!((distance-0.0001)>FRET_cutoff)){
                    hops_temp[index].setObjectIt(object_it);
                    hops_temp[index].setDestCoords(dest_coords);
                    E_delta = (getSiteEnergy(dest_coords)-getSiteEnergy(object_coords));
                    hops_temp[index].calculateExecutionTime(R_exciton_hopping,distance,E_delta,getTemperature(),getTime());
                    hops_valid[index] = true;
                }
                else{
                    hops_valid[index] = false;
                }
            }
        }
    }
    // Exciton Recombination
    auto recombination_event_it = exciton_recombination_events.begin();
    advance(recombination_event_it,std::distance(excitons.begin(),exciton_it));
    rate = 1/Exciton_lifetime;
    recombination_event_it->calculateExecutionTime(rate,getTime());
    // Determine the fastest available hop event
    bool No_hops_valid = true;
    auto hop_target_it = hops_temp.end();
    for(auto hop_it=hops_temp.begin();hop_it!=hops_temp.end();++hop_it){
        if(hops_valid[std::distance(hops_temp.begin(),hop_it)] && (hop_target_it==hops_temp.end() || hop_it->getExecutionTime()<hop_target_it->getExecutionTime())){
            hop_target_it = hop_it;
            No_hops_valid = false;
        }
    }
    // Compare fastest hop event with recombination event to determine fastest event for this exciton
    if(!No_hops_valid && hop_target_it->getExecutionTime() < recombination_event_it->getExecutionTime()){
        auto hop_list_it = exciton_hop_events.begin();
        advance(hop_list_it,std::distance(excitons.begin(),exciton_it));
        *hop_list_it = *hop_target_it;
        setEvent(exciton_it->getEventIt(),&(*hop_list_it));
    }
    else{
        setEvent(exciton_it->getEventIt(),&(*recombination_event_it));
    }
}

void Exciton_sim::calculateExcitonListEvents(const vector<list<Object*>::iterator>& exciton_it_vec){
    for(auto it=exciton_it_vec.begin();it!=exciton_it_vec.end();++it){
        calculateExcitonEvents(*it);
    }
}

bool Exciton_sim::checkFinished(){
    if(Enable_diffusion_test){
        return (N_excitons_recombined==N_tests);
    }
    cout << "Error checking simulation finish conditions." << endl;
    return true;
}

void Exciton_sim::deleteExciton(const list<Object*>::iterator object_it){
    auto exciton_it = getExcitonIt(*object_it);
    // remove the Object and Event pointers from Simulation
    removeObject(object_it);
    // Locate corresponding recombine event
    auto recombination_list_it = exciton_recombination_events.begin();
    advance(recombination_list_it,std::distance(excitons.begin(),exciton_it));
    // Locate corresponding hop event
    auto hop_list_it = exciton_hop_events.begin();
    advance(hop_list_it,std::distance(excitons.begin(),exciton_it));
    // Delete exciton
    excitons.erase(exciton_it);
    // Delete exciton recombination event
    exciton_recombination_events.erase(recombination_list_it);
    // Delete exciton hop event
    exciton_hop_events.erase(hop_list_it);
}

bool Exciton_sim::executeExcitonCreation(const list<Event*>::iterator event_it){
    // Determine coords for the new exciton
    const Coords coords_new = calculateExcitonCreationCoords();
    // Create the new exciton and at it to the simulation
    Exciton exciton_new(getTime(),N_excitons_created+1,coords_new);
    excitons.push_back(exciton_new);
    auto object_it = addObject(&excitons.back());
    // Add an empty hop and recombine event to the corresponding lists
    Exciton_Hop hop_event;
    exciton_hop_events.push_back(hop_event);
    Exciton_Recombination recombination_event;
    recombination_event.setObjectIt(object_it);
    exciton_recombination_events.push_back(recombination_event);
    // Update exciton counters
    N_excitons_created++;
    N_excitons++;
    // Log event
    if(loggingEnabled()){
        *Logfile << "Created exciton " << exciton_new.getTag() << " at site " << coords_new.x << "," << coords_new.y << "," << coords_new.z << "." << endl;
    }
    // Find all nearby excitons and calculate their events
    auto neighbors = findRecalcNeighbors(coords_new);
    calculateExcitonListEvents(neighbors);
    // Calculate next exciton creation event
    exciton_creation_event.calculateExecutionTime(R_exciton_generation,getTime());
    return true;
}

bool Exciton_sim::executeExcitonHop(const list<Event*>::iterator event_it){
    if(isOccupied((*event_it)->getDestCoords())){
        cout << "Exciton hop cannot be executed. Destination site is already occupied." << endl;
        return false;
    }
    else{
        // Get event info
        Coords coords_initial = (*((*event_it)->getObjectIt()))->getCoords();
        // Move the exciton in the Simulation
        moveObject((*event_it)->getObjectIt(),(*event_it)->getDestCoords());
        // Log event
        if(loggingEnabled()){
            *Logfile << "Exciton " << (*((*event_it)->getObjectIt()))->getTag() << " hopped to site " << (*event_it)->getDestCoords().x << "," << (*event_it)->getDestCoords().y << "," << (*event_it)->getDestCoords().z << "." << endl;
        }
        // Find all nearby excitons and calculate their events
        vector<list<Object*>::iterator> neighbors = findRecalcNeighbors(coords_initial);
        vector<list<Object*>::iterator> neighbors2 = findRecalcNeighbors((*event_it)->getDestCoords());
        neighbors.insert(neighbors.end(),neighbors2.begin(),neighbors2.end());
        removeObjectItDuplicates(neighbors);
        calculateExcitonListEvents(neighbors);
        return true;
    }
}

bool Exciton_sim::executeExcitonRecombination(const list<Event*>::iterator event_it){
    // Get event info
    int exciton_tag = (*((*event_it)->getObjectIt()))->getTag();
    Coords coords_initial = (*((*event_it)->getObjectIt()))->getCoords();
    // Output diffusion distance
    if(Enable_diffusion_test){
        diffusion_distances.push_back(getUnitSize()*(*((*event_it)->getObjectIt()))->calculateDisplacement());
    }
    // delete Exciton and its events in the Exciton_sim class
    deleteExciton((*event_it)->getObjectIt());
    // Update exciton counters
    N_excitons--;
    N_excitons_recombined++;
    // Log event
    if(loggingEnabled()){
        *Logfile << "Exciton " << exciton_tag << " recombined at site " << coords_initial.x << "," << coords_initial.y << "," << coords_initial.z << "." << endl;
    }
    // Find all nearby excitons and calculate their events
    vector<list<Object*>::iterator> neighbors = findRecalcNeighbors(coords_initial);
    calculateExcitonListEvents(neighbors);
    return true;
}

bool Exciton_sim::executeNextEvent(){
    auto event_it = chooseNextEvent();
    string event_name = (*event_it)->getName();
    if(loggingEnabled()){
        *Logfile << "Executing " << event_name << " event" << endl;
    }
    // Update simulation time
    updateTime((*event_it)->getExecutionTime());
    // Determine which event has been chosen
    if(event_name.compare(Exciton_Creation::name)==0){
        return executeExcitonCreation(event_it);
    }
    else if(event_name.compare(Exciton_Hop::name)==0){
        return executeExcitonHop(event_it);
    }
    else if(event_name.compare(Exciton_Recombination::name)==0){
        return executeExcitonRecombination(event_it);
    }
    else{
        //error
        cout << "Valid event not found when calling executeNextEvent" << endl;
        return false;
    }
}

list<Exciton>::iterator Exciton_sim::getExcitonIt(const Object* object_ptr){
    for(auto exciton_it=excitons.begin();exciton_it!=excitons.end();++exciton_it){
        if(object_ptr->getTag()==exciton_it->getTag()){
            return exciton_it;
        }
    }
}

vector<double> Exciton_sim::getDiffusionData(){
    return diffusion_distances;
}

int Exciton_sim::getN_excitons_created(){
    return N_excitons_created;
}

int Exciton_sim::getN_excitons_recombined(){
    return N_excitons_recombined;
}

void Exciton_sim::outputStatus(){
    cout << getId() << ": Time = " << getTime() << " seconds.\n";
    cout << getId() << ": " << N_excitons_created << " excitons have been created and " << getN_events_executed() << " events have been executed.\n";
    cout << getId() << ": There are " << N_excitons << " excitons in the lattice.\n";
    for(auto it=excitons.begin();it!=excitons.end();++it){
        cout << getId() << ": Exciton " << it->getTag() << " is at " << it->getCoords().x << "," << it->getCoords().y << "," << it->getCoords().z << ".\n";
    }
    cout.flush();
}

double Exciton_sim::getSiteEnergy(const Coords& coords){
    return sites[getSiteIndex(coords)].getEnergy();
}
