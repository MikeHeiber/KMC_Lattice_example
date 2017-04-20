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
    for(int i=0;i<lattice.size();i++){
        if(Enable_gaussian_dos || Enable_exponential_dos){
            site.setEnergyIt(site_energies.begin()+i);
        }
        else{
            site.setEnergyIt(site_energies.begin());
        }
        sites.push_back(site);
        lattice.push_back(unique_ptr<Site>(&sites.back()));
    }
    // Initialize exciton creation event
    calculateExcitonCreationEvent();
    events.push_back(unique_ptr<Event>(&exciton_creation_event));
    exciton_creation_it = --events.end();
}

double Exciton_sim::calculateDiffusionLength_avg(){
    return vector_avg(diffusion_distances);
}

double Exciton_sim::calculateDiffusionLength_stdev(){
    return vector_stdev(diffusion_distances);
}

void Exciton_sim::calculateExcitonCreationEvent(){
    exciton_creation_event.calculateEvent(Exciton_generation_rate);
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
                    event_hop.calculateEvent(dest_coords,distance,E_delta,getTemperature());
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
    event_recombine.calculateEvent();
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

bool Exciton_sim::executeNextEvent(){
    (*chooseNextEvent())->executeEvent();

}

int Exciton_sim::getN_excitons_created(){
    return N_excitons_created;
}

double Exciton_sim::getSiteEnergy(const Coords& coords){
    auto derived_ptr = static_unique_ptr_cast<Site_OSC>(std::move(*getSiteIt(coords)));
    return derived_ptr->getEnergy();
}
