// Copyright (c) 2017 Michael C. Heiber
// This source file is part of the KMC_Lattice_example project, which is subject to the MIT License.
// For more information, see the LICENSE file that accompanies this software.
// The KMC_Lattice_example project can be found on Github at https://github.com/MikeHeiber/KMC_Lattice_example

#ifndef EXCITON_SIM_H
#define EXCITON_SIM_H

#include "KMC_Lattice/Utils.h"
#include "KMC_Lattice/Simulation.h"
#include "KMC_Lattice/Object.h"
#include "KMC_Lattice/Event.h"
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
        void setEnergyIt(const vector<double>::iterator it){energy_it = it;}
        double getEnergy(){return *energy_it;}
    private:
        vector<double>::iterator energy_it;
};

class Exciton_sim : public Simulation{
    public:
        // Functions
        Exciton_sim(const Parameters_Exciton_Sim& params,const int id);
        double calculateDiffusionLength_avg();
        double calculateDiffusionLength_stdev();
        bool checkFinished();
        bool executeNextEvent();
        vector<double> getDiffusionData();
        int getN_excitons_created();
        int getN_excitons_recombined();
        void outputStatus();
    protected:
        vector<Site_OSC> sites;
        list<Exciton> excitons;
        list<Exciton_Hop> exciton_hop_events;
        list<Exciton_Recombination> exciton_recombination_events;
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
        // Derived Paramters
        double R_exciton_generation;
        // Additional Data Structures
        vector<double> site_energies;
        vector<double> diffusion_distances;
        Exciton_Creation exciton_creation_event;
        list<Event*>::iterator exciton_creation_it;
        // Additional Counters
        int N_excitons_created;
        int N_excitons_recombined;
        int N_excitons;
        // Additional Functions
        Coords calculateExcitonCreationCoords();
        void calculateExcitonEvents(const list<Object*>::iterator object_it);
        void calculateExcitonListEvents(const vector<list<Object*>::iterator>& exciton_it_vec);
        void deleteExciton(const list<Object*>::iterator object_it);
        bool executeExcitonCreation(const list<Event*>::iterator event_it);
        bool executeExcitonHop(const list<Event*>::iterator event_it);
        bool executeExcitonRecombination(const list<Event*>::iterator event_it);
        list<Exciton>::iterator getExcitonIt(const Object* object_ptr);
        double getSiteEnergy(const Coords& coords);
};

#endif // EXCITON_SIM_H
