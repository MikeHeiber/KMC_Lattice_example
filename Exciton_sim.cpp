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
    // Initialize counters
    N_excitons = 0;
    N_excitons_created = 0;
    // Initialize sites
    sites.clear();
    Site_OSC site;
    site.clearOccupancy();
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
    double R_exciton_generation = Exciton_generation_rate*getNumSites()*intpow(1e-7*getUnitSize(),3);
    exciton_creation_event.calculateEvent(dest_coords,0,0,getTemperature(),R_exciton_generation);
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
        cout << "Attempting to create exciton at " << dest_coords.x << "," << dest_coords.y << "," << dest_coords.z << "." << endl;
        if(loggingEnabled()){
            ostringstream msg;
            msg << "Attempting to create exciton at " << dest_coords.x << "," << dest_coords.y << "," << dest_coords.z << "." << endl;
            logMSG(msg);
        }
        if(isOccupied(dest_coords)){
            cout << "Site is already occupied." << endl;
            continue;
        }
        else{
            success = true;
            return dest_coords;
        }
    }
}

void Exciton_sim::calculateExcitonEvents(const list<unique_ptr<Object>>::iterator object_it){
    cout << "Calculating events for exciton " << (*object_it)->getTag() << "." << endl;
    const list<Exciton>::iterator exciton_it = getExcitonIt(*object_it);
    Coords object_coords = exciton_it->getCoords();
    Coords dest_coords;
    double distance,E_delta;
    int index;
    // Exciton Hopping
    Exciton_Hop event_hop;
    event_hop.setObjectIt(object_it);
    static const int range = ceil(FRET_cutoff/getUnitSize());
    static const int dim = (2*range+1);
    static vector<Exciton_Hop> hops_temp(dim*dim*dim,event_hop);
    static vector<bool> hops_valid(dim*dim*dim,false);
    for(int i=-range;i<=range;i++){
        for(int j=-range;j<=range;j++){
            for(int k=-range;k<=range;k++){
                if(i==0 && j==0 && k==0){
                    hops_valid[(i+range)*dim*dim+(j+range)*dim+(k+range)] = false;
                    continue;
                }
                dest_coords.x = object_coords.x+i+calculateDX(object_coords.x,i);
                dest_coords.y = object_coords.y+j+calculateDY(object_coords.y,j);
                dest_coords.z = object_coords.z+k+calculateDZ(object_coords.z,k);
                distance = getUnitSize()*sqrt((double)(i*i+j*j+k*k));
                index = (i+range)*dim*dim+(j+range)*dim+(k+range);
                if(!((distance-0.0001)>FRET_cutoff) && !((*getSiteIt(dest_coords))->isOccupied())){
                    hops_temp[index].setObjectIt(object_it);
                    E_delta = (getSiteEnergy(dest_coords)-getSiteEnergy(object_coords));
                    hops_temp[index].calculateEvent(dest_coords,distance,E_delta,getTemperature(),0);
                    hops_valid[index] = true;
                }
                else{
                    hops_valid[index] = false;
                }
            }
        }
    }
    // Exciton Recombination
    Exciton_Recombine event_recombine;
    event_recombine.setObjectIt(object_it);
    event_recombine.calculateEvent(object_coords,0,0,0,0);
    list<Exciton_Recombine>::iterator recombine_list_it = exciton_recombine_events.begin();
    advance(recombine_list_it,std::distance(excitons.begin(),exciton_it));
    *recombine_list_it = event_recombine;
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
        cout << "Hop event was faster than the recombination event." << endl;
        hop_list_it = exciton_hop_events.begin();
        advance(hop_list_it,std::distance(excitons.begin(),exciton_it));
        *hop_list_it = *hop_target_it;
        setEvent(exciton_it->getEventIt(),unique_ptr<Event>(&(*hop_list_it)));
    }
    else{
        cout << "Recombination event was faster than the hop event." << endl;
        setEvent(exciton_it->getEventIt(),unique_ptr<Event>(&(*recombine_list_it)));
    }
}

bool Exciton_sim::checkFinished(){
    if(Enable_diffusion_test){
        return (N_excitons_created>N_tests);
    }
}

void Exciton_sim::deleteExciton(const list<Exciton>::iterator exciton_it){
    // locate corresponding recombine event
    list<Exciton_Recombine>::iterator recombine_list_it = exciton_recombine_events.begin();
    advance(recombine_list_it,std::distance(excitons.begin(),exciton_it));
    // locate corresponding hop event
    list<Exciton_Hop>::iterator hop_list_it = exciton_hop_events.begin();
    advance(hop_list_it,std::distance(excitons.begin(),exciton_it));
    // delete exciton
    excitons.erase(exciton_it);
    // delete exciton recombination event
    exciton_recombine_events.erase(recombine_list_it);
    // delete exciton hop event
    exciton_hop_events.erase(hop_list_it);
}

bool Exciton_sim::executeExcitonCreation(const list<unique_ptr<Event>>::iterator event_it){
    incrementTime((*event_it)->getWaitTime());
    Coords coords_new = calculateExcitonCreationCoords();
    Exciton exciton_new(getTime(),N_excitons_created+1,coords_new);
    excitons.push_back(exciton_new);
    unique_ptr<Object> object_ptr = unique_ptr<Object>(&excitons.back());
    list<unique_ptr<Object>>::iterator object_it = addObject(object_ptr);
    cout << "Created exciton " << exciton_new.getTag() << " at " << coords_new.x << "," << coords_new.y << "," << coords_new.z << " at " << getTime() << " seconds." << endl;
    N_excitons_created++;
    N_excitons++;
    Exciton_Hop hop_event;
    exciton_hop_events.push_back(hop_event);
    Exciton_Recombine recombine_event;
    exciton_recombine_events.push_back(recombine_event);
    calculateExcitonEvents(object_it);
    if(loggingEnabled()){
        ostringstream msg;
        msg << "Created exciton " << exciton_new.getTag() << " at " << coords_new.x << "," << coords_new.y << "," << coords_new.z << " at " << getTime() << " seconds." << endl;
        logMSG(msg);
    }
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
        calculateExcitonEvents((*event_it)->getObjectIt());
        if(loggingEnabled()){
            ostringstream msg;
            msg << "Exciton " << (*((*event_it)->getObjectIt()))->getTag() << " hopped to " << (*event_it)->getDestCoords().x << "," << (*event_it)->getDestCoords().y << "," << (*event_it)->getDestCoords().z << "." << endl;
            logMSG(msg);
        }
        return true;
    }
}

bool Exciton_sim::executeExcitonRecombine(const list<unique_ptr<Event>>::iterator event_it){
    incrementTime((*event_it)->getWaitTime());
    deleteExciton(getExcitonIt(*((*event_it)->getObjectIt())));
    removeObject((*event_it)->getObjectIt());
    N_excitons--;
    return true;
}

bool Exciton_sim::executeNextEvent(){
    list<unique_ptr<Event>>::iterator event_it = chooseNextEvent();
    string event_name = (*event_it)->getName();
    cout << event_name << " event chosen with wait time = " << (*event_it)->getWaitTime() << " s." << endl;
    if(loggingEnabled()){
        ostringstream msg;
        msg << "Executing " << event_name << " event" << endl;
        logMSG(msg);
    }
    if(event_name.compare(Exciton_Creation::name)==0){
        return executeExcitonCreation(event_it);
    }
    else if(event_name.compare(Exciton_Hop::name)==0){
        return executeExcitonHop(event_it);
    }
    else if(event_name.compare(Exciton_Recombine::name)==0){
        return executeExcitonRecombine(event_it);
    }
    else{
        //error
        cout << "Valid event not found when calling executeNextEvent" << endl;
        return false;
    }
}

list<Exciton>::iterator Exciton_sim::getExcitonIt(const unique_ptr<Object>& object_ptr){
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

float Exciton_sim::getSiteEnergy(const Coords& coords){
    return sites[getSiteIndex(coords)].getEnergy();
}
