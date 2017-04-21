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
    Exciton::R_hop = R_exciton_hopping;
    Exciton::lifetime = Exciton_lifetime;
    Exciton::FRET_cutoff = FRET_cutoff;
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
        createGaussianDOSVector(site_energies,0,Site_energy_stdev,getId());
    }
    else if(Enable_exponential_dos){
        site_energies.assign(getNumSites(),0);
        createExponentialDOSVector(site_energies,0,Site_energy_urbach,getId());
    }
    else{
        site_energies.push_back(0);
    }
    // Initialize sites
    sites.clear();
    Site_OSC site;
    for(int i=0;i<getNumSites();i++){
        if(Enable_gaussian_dos || Enable_exponential_dos){
            site.setEnergyIt(site_energies.begin()+i);
        }
        else{
            site.setEnergyIt(site_energies.begin());
        }
        sites.push_back(site);
        unique_ptr<Site> site_ptr = unique_ptr<Site>(&sites.back());
        addSite(site_ptr);
    }
    // Initialize exciton creation event
    Coords dest_coords;
    exciton_creation_event.calculateEvent(dest_coords,0,0,getTemperature(),Exciton_generation_rate);
    unique_ptr<Event> event_ptr = unique_ptr<Event>(&exciton_creation_event);
    exciton_creation_it = addEvent(event_ptr);
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
            stringstream msg;
            msg << "Attempting to create exciton at " << dest_coords.x << "," << dest_coords.y << "," << dest_coords.z << "." << endl;
            logMSG(msg);
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

void Exciton_sim::calculateExcitonEvents(const list<Exciton>::iterator exciton_it){
    Coords object_coords = exciton_it->getCoords();
    Coords dest_coords;
    double distance,E_delta;
    // Exciton Hopping
    Exciton_Hop event_hop;
    static const int range = ceil(FRET_cutoff/getUnitSize());
    static vector<Exciton_Hop> hops_temp(range*range*range);
    static vector<bool> hops_valid(range*range*range);
    for(int i=-range;i<=range;i++){
        for(int j=-range;j<=range;j++){
            for(int k=-range;k<=range;k++){
                dest_coords.x = object_coords.x+i+calculateDX(object_coords.x,i);
                dest_coords.y = object_coords.y+j+calculateDY(object_coords.y,j);
                dest_coords.z = object_coords.z+k+calculateDZ(object_coords.z,k);
                distance = getUnitSize()*sqrt((double)(i*i+j*j+k*k));
                if(!(*getSiteIt(dest_coords))->isOccupied() && !distance>FRET_cutoff){
                    E_delta = getSiteEnergy(dest_coords)-getSiteEnergy(object_coords);
                    event_hop.calculateEvent(dest_coords,distance,E_delta,getTemperature(),0);
                    hops_temp[i*range*range+j*range+k] = event_hop;
                    hops_valid[i*range*range+j*range+k] = true;
                }
                else{
                    hops_valid[i*range*range+j*range+k] = false;
                }
            }
        }
    }
    // Exciton Recombination
    Exciton_Recombine event_recombine;
    event_recombine.calculateEvent(object_coords,0,0,0,0);
    // Determine the fastest available hop event
    vector<Exciton_Hop>::iterator hop_it;
    vector<Exciton_Hop>::iterator hop_target_it = hops_temp.end();
    for(hop_it=hops_temp.begin();hop_it!=hops_temp.end();++hop_it){
        if(hops_valid[std::distance(hops_temp.begin(),hop_it)] && (hop_target_it==hops_temp.end() || hop_it->getWaitTime()<hop_target_it->getWaitTime())){
            hop_target_it = hop_it;
        }
    }
    // Compare fastest hop event with recombination event to determine fastest event for this exciton
    list<Exciton_Hop>::iterator hop_list_it;
    if(hop_target_it->getWaitTime() < event_recombine.getWaitTime()){
        hop_list_it = exciton_hop_events.begin();
        advance(hop_list_it,std::distance(excitons.begin(),exciton_it));
        *hop_list_it = *hop_target_it;
        setEvent(exciton_it->getEventIt(),unique_ptr<Event>(&(*hop_list_it)));
    }
    else{
        setEvent(exciton_it->getEventIt(),unique_ptr<Event>(&event_recombine));
    }
}

bool Exciton_sim::checkFinished(){
    if(Enable_diffusion_test){
        return (N_excitons_created>N_tests);
    }
}

bool Exciton_sim::executeExcitonCreation(const list<unique_ptr<Event>>::iterator event_it){
    incrementTime((*event_it)->getWaitTime());
    Coords coords_new = calculateExcitonCreationCoords();
    Exciton exciton_new(getTime(),N_excitons_created+1,coords_new);
    excitons.push_back(exciton_new);
    unique_ptr<Object> object_ptr = unique_ptr<Object>(&excitons.back());
    addObject(object_ptr);
    N_excitons_created++;
    N_excitons++;
    calculateExcitonEvents(--excitons.end());
    return true;
}

bool Exciton_sim::executeExcitonHop(const list<unique_ptr<Event>>::iterator event_it){
    if(isOccupied((*event_it)->getDestCoords())){
        cout << "Exciton hop cannot be executed. Destination site is already occupied." << endl;
        return false;
    }
    else{
        incrementTime((*event_it)->getWaitTime());
        moveObject((*event_it)->getObjectIt(),(*event_it)->getDestCoords());
        calculateExcitonEvents(getExcitonIt(*((*event_it)->getObjectIt())));
        return true;
    }
}

bool Exciton_sim::executeExcitonRecombine(const list<unique_ptr<Event>>::iterator event_it){
    incrementTime((*event_it)->getWaitTime());
    deleteExciton((*event_it)->getObjectIt());
    removeObject((*event_it)->getObjectIt());
    N_excitons--;
    return true;
}

bool Exciton_sim::executeNextEvent(){
    list<unique_ptr<Event>>::iterator event_it = chooseNextEvent();
    string event_name = (*event_it)->getName();
    if(event_name.compare("Exciton Creation")){
        return executeExcitonCreation(event_it);
    }
    else if(event_name.compare("Exciton Hop")){
        return executeExcitonHop(event_it);
    }
    else if(event_name.compare("Exciton Recombine")){
        return executeExcitonRecombine(event_it);
    }
    else{
        //error
        return false;
    }
}

list<Exciton>::iterator Exciton_sim::getExcitonIt(unique_ptr<Object>& object_ptr){
    list<Exciton>::iterator exciton_it;
    for(exciton_it=excitons.begin();exciton_it!=excitons.end();++exciton_it){
        if(object_ptr->getTag()==exciton_it->getTag()){
            return exciton_it;
        }
    }
}

int Exciton_sim::getN_excitons_created(){
    return N_excitons_created;
}

double Exciton_sim::getSiteEnergy(const Coords& coords){
    auto derived_ptr = static_unique_ptr_cast<Site_OSC>(std::move(*getSiteIt(coords)));
    return derived_ptr->getEnergy();
}
