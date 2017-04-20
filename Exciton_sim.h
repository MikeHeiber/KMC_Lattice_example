#ifndef EXCITON_SIM_H
#define EXCITON_SIM_H

#include "Utils.h"
#include "Simulation.h"
#include "Object.h"
#include "Event.h"
#include "Exciton.h"

using namespace std;

struct Parameters_Exciton_Sim : Parameters_Simulation{
    // General Parameters
    // Test Parameters
    bool Enable_diffusion_test;
    int N_tests;
    // Object Parameters
    // Excitons
    double Exciton_generation_rate; // 1/seconds
    double Exciton_lifetime; // seconds
    double R_exciton_hopping;
    int FRET_cutoff; // nm
    // Lattice Site Parameters
    bool Enable_gaussian_dos;
    double Site_energy_stdev; // eV
    bool Enable_exponential_dos;
    double Site_energy_urbach;
    // Output files
};

class Site_OSC : public Site{
    public:
        void setEnergyIt(const vector<float>::iterator it){energy_it = it;}
        double getEnergy(){return *energy_it;}
    private:
        vector<float>::iterator energy_it;
};

class Exciton_sim : public Simulation{
    public:
        // Functions
        Exciton_sim(const Parameters_Exciton_Sim& params,const int id);
        double calculateDiffusionLength_avg();
        double calculateDiffusionLength_stdev();
        bool checkFinished();
        bool executeNextEvent();
        int getN_excitons_created();
        double getSiteEnergy(const Coords& coords);
    protected:
        vector<Site_OSC> sites;
        list<Exciton> excitons;
        list<Exciton_Hop> exciton_hop_events;
        list<Exciton_Recombine> exciton_recombine_events;
    private:
        // Additional General Parameters
        //
        // Test Parameters
        bool Enable_diffusion_test;
        int N_tests;
        // Object Parameters
        double Exciton_generation_rate;
        double Exciton_lifetime; // seconds
        double R_exciton_hopping;
        int FRET_cutoff;
        // Additional Lattice Parameters
        bool Enable_gaussian_dos;
        double Site_energy_stdev; // eV
        bool Enable_exponential_dos;
        double Site_energy_urbach;
        // Additional Output Files
        //
        // Additional Data Structures
        vector<float> site_energies;
        vector<double> diffusion_distances;
        Exciton_Creation exciton_creation_event;
        list<unique_ptr<Event>>::iterator exciton_creation_it;
        // Additional Counters
        int N_excitons_created;
        int N_excitons;
        // Additional Functions
        void calculateExcitonCreationEvent();
        void calculateExcitonEvents(const list<Exciton>::iterator object_it);
};

#endif // EXCITON_SIM_H
