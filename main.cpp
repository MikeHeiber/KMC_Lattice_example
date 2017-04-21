#include "Utils.h"
#include "Exciton_sim.h"
#include "mpi.h"
#include <fstream>
#include <string>
#include <vector>
//#include <iostream>
//#include <sstream>
//#include <cmath>

using namespace std;

struct Parameters_main{
    bool Enable_mpi;
};

//Declare Functions
bool importParameters(ifstream * inputfile,Parameters_main& params_main,Parameters_Exciton_Sim& params);

int main(int argc,char *argv[]){
    string version = "v0.1-alpha";
    // Parameters
    bool End_sim = false;
    // File declaration
    ifstream parameterfile;
    ofstream resultsfile;
    ofstream logfile;
    ofstream outputfile;
    stringstream ss;
    // Initialize variables
    string parameterfilename;
    string logfilename;
    Parameters_main params_main;
    Parameters_Exciton_Sim params_sim;
    int nproc = 1;
    int procid = 0;
    int elapsedtime;
    time_t start, end;
    bool eventResult;
    bool success;
    // Start timer
    start = time(NULL);
    // Import parameters and options from file and command line arguments
    cout << "Loading input parameters from file... ";
    parameterfilename = argv[1];
    parameterfile.open(parameterfilename.c_str(),ifstream::in);
    if(!parameterfile){
        cout << "Error loading parameter file.  Program will now exit." << endl;
        return 0;
    }
    success = importParameters(&parameterfile,params_main,params_sim);
    parameterfile.close();
    if(!success){
        cout << "Error importing parameters from parameter file.  Program will now exit." << endl;
        return 0;
    }
    cout << "Parameter loading complete!" << endl;
    // Initialize mpi options
    if(params_main.Enable_mpi){
        cout << "Initializing MPI options... ";
        MPI::Init(argc,argv);
        nproc = MPI::COMM_WORLD.Get_size();
        procid = MPI::COMM_WORLD.Get_rank();
        cout << "MPI initialization complete!" << endl;
    }
    // Setup file output
    cout << procid << ": Creating output files..." << endl;
    if(params_sim.Enable_logging){
        ss << "log" << procid << ".txt";
        logfilename = ss.str();
        logfile.open(ss.str().c_str());
        ss.str("");
    }
    ss << "output" << procid << ".txt";
    outputfile.open(ss.str().c_str());
    ss.str("");
    // Add file input and output information to input parameters
    params_sim.Logfile = &logfile;
    // Initialize Simulation
    cout << procid << ": Initializing simulation " << procid << "..." << endl;
    Exciton_sim sim(params_sim,procid);
    cout << procid << ": Simulation initialization complete" << endl;
    // Output initial status
    if(params_sim.Enable_logging){
        //sim.outputStatus(&logfile);
    }
    // Begin Simulation loop
    while(!End_sim){
        eventResult = sim.executeNextEvent();
        if(!eventResult){
            cout << procid << ": Event execution failed, simulation will now terminate." << endl;
            break;
        }
        // Check if simulation has finished
        End_sim = sim.checkFinished();
        // Output status
        if(params_sim.Enable_logging){
            //sim.outputStatus(&logfile);
            //sim.outputEventQueue(5,&logfile);
            if(sim.getN_events_executed()%2000==0){
                //logfile.close();
                //logfile.open(logfilename.c_str());
            }
        }
    }
    // Calculate final statistics
    cout << procid << ": Simulation finished." << endl;
    end = time(NULL);
    elapsedtime = difftime(end,start);
    // Output calculation parameters
    double exciton_diffusion_avg = 0;
    double exciton_diffusion_stdev = 0;
    if(params_sim.Enable_diffusion_test){
        exciton_diffusion_avg = sim.calculateDiffusionLength_avg();
        exciton_diffusion_stdev = sim.calculateDiffusionLength_stdev();
    }
    ss << "results" << procid << ".txt";
    resultsfile.open(ss.str().c_str());
    ss.str("");
    resultsfile << "KMC_Lattice " << version << ":\n";
    resultsfile << "Simulation on processor " << procid << " finished.\n";
    resultsfile << "Calculation time elapsed is " << (double)elapsedtime/60 << " minutes.\n";
    resultsfile << sim.getTime() << " seconds have been simulated.\n";
    resultsfile << sim.getN_events_executed() << " events have been executed.\n";
    resultsfile << sim.getN_excitons_created() << " objects have been created.\n";
    if(params_sim.Enable_diffusion_test){
        resultsfile << "Exciton diffusion test results:\n";
        resultsfile << "Effective Exciton Diffusion Length is " << exciton_diffusion_avg << " ± " << exciton_diffusion_stdev << " nm\n";
    }
    resultsfile << endl;
//    if(Enable_mpi){
//        if(params.Enable_diffusion_test){
//            // Generate Results from Polaron log
//            if(procid==0){
//                cout << "Analyzing diffusion log files." << endl;
//            }
//            int N = (12+log10(params.Time_cutoff))*20;
//            int array_size = N+1;
//            double* exciton_displacement_data;
//            double* exciton_energy_data;
//            int* exciton_count_data;
//            double* time_data = (double*)malloc(sizeof(double)*array_size);
//            exciton_displacement_data = (double*)malloc(sizeof(double)*array_size);
//            exciton_energy_data = (double*)malloc(sizeof(double)*array_size);
//            exciton_count_data = (int*)malloc(sizeof(int)*array_size);
//            vector<double> times_temp;
//            vector<double> displacements_temp;
//            vector<double> energies_temp;
//            string line,var;
//            int n,tag,type,tag_target,type_target;
//            double time,displacement,energy;
//            bool first_loop;
//            // Load diffusion log file for one proc at a time
//            for(int i=0;i<nproc;i++){
//                MPI_Barrier(MPI_COMM_WORLD);
//                if(i==procid){
//                    first_loop = true;
//                    for(int i=0;i<=N;i++){
//                        time_data[i] = pow(10,(0.05*i-12));
//                        exciton_displacement_data[i] = 0;
//                        exciton_energy_data[i] = 0;
//                        exciton_count_data[i] = 0;
//                    }
//                }
//            }
//            double* exciton_displacement_data_all;
//            double* exciton_energy_data_all;
//            int *exciton_count_data_all;
//            if(procid==0){
//                exciton_displacement_data_all = (double *)malloc(sizeof(double)*nproc*array_size);
//                exciton_energy_data_all = (double *)malloc(sizeof(double)*nproc*array_size);
//                exciton_count_data_all = (int *)malloc(sizeof(int)*nproc*array_size);
//            }
//            MPI_Barrier(MPI_COMM_WORLD);
//            MPI_Gather(exciton_displacement_data,array_size,MPI_DOUBLE,exciton_displacement_data_all,array_size,MPI_DOUBLE,0,MPI_COMM_WORLD);
//            MPI_Gather(exciton_energy_data,array_size,MPI_DOUBLE,exciton_energy_data_all,array_size,MPI_DOUBLE,0,MPI_COMM_WORLD);
//            MPI_Gather(exciton_count_data,array_size,MPI_INT,exciton_count_data_all,array_size,MPI_INT,0,MPI_COMM_WORLD);
//            if(procid==0){
//                ofstream displacementfile;
//                displacementfile.open("displacementlog.txt");
//                double exciton_displacement_avg,exciton_energy_avg;
//                int exciton_count;
//                for(int i=0;i<array_size;i++){
//                    exciton_displacement_avg = 0;
//                    exciton_energy_avg = 0;
//                    exciton_count = 0;
//                    for(int j=0;j<nproc;j++){
//                        exciton_displacement_avg += exciton_displacement_data_all[j*array_size+i];
//                        exciton_energy_avg += exciton_energy_data_all[j*array_size+i];
//                        exciton_count += exciton_count_data_all[j*array_size+i];
//                    }
//                    displacementfile << time_data[i] << "," << exciton_displacement_avg/exciton_count << "," << exciton_energy_avg/exciton_count << endl;
//                }
//                displacementfile.close();
//            }
//        }
//    }
    // Cleanup
    cout << procid << ": Begin Cleanup" << endl;
    if(params_sim.Enable_logging){
        logfile.close();
    }
    resultsfile.close();
    outputfile.close();
    if(params_main.Enable_mpi){
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Finalize();
    }
    return 0;
}

bool importParameters(ifstream * inputfile,Parameters_main& params_main,Parameters_Exciton_Sim& params){
    string line;
    string var;
    size_t pos;
    vector<string> stringvars;
    while(inputfile->good()){
        getline(*inputfile,line);
        if((line.substr(0,2)).compare("--")!=0 && (line.substr(0,2)).compare("##")!=0){
            pos = line.find("/",0);
            var = line.substr(0,pos-1);
            stringvars.push_back(var);
        }
    }
    int i = 0;
    // General Parameters
    //enable_mpi
    if(stringvars[i].compare("true")==0){
        params_main.Enable_mpi = true;
    }
    else if(stringvars[i].compare("false")==0){
        params_main.Enable_mpi = false;
    }
    else{
        cout << "Error setting mpi options" << endl;
        return false;
    }
    i++;
    //enable_logging
    if(stringvars[i].compare("true")==0){
        params.Enable_logging = true;
    }
    else if(stringvars[i].compare("false")==0){
        params.Enable_logging = false;
    }
    else{
        cout << "Error setting logging options" << endl;
        return false;
    }
    i++;
    //enable_periodic_x
    if(stringvars[i].compare("true")==0){
        params.Enable_periodic_x = true;
    }
    else if(stringvars[i].compare("false")==0){
        params.Enable_periodic_x = false;
    }
    else{
        cout << "Error setting x-periodic boundary options" << endl;
        return false;
    }
    i++;
    //enable_periodic_y
    if(stringvars[i].compare("true")==0){
        params.Enable_periodic_y = true;
    }
    else if(stringvars[i].compare("false")==0){
        params.Enable_periodic_y = false;
    }
    else{
        cout << "Error setting y-periodic boundary options" << endl;
        return false;
    }
    i++;
    //enable_periodic_z
    if(stringvars[i].compare("true")==0){
        params.Enable_periodic_z = true;
    }
    else if(stringvars[i].compare("false")==0){
        params.Enable_periodic_z = false;
    }
    else{
        cout << "Error setting z-periodic boundary options" << endl;
        return false;
    }
    i++;
    params.Length = atoi(stringvars[i].c_str());
    i++;
    params.Width = atoi(stringvars[i].c_str());
    i++;
    params.Height = atoi(stringvars[i].c_str());
    i++;
    params.Unit_size = atof(stringvars[i].c_str());
    i++;
    params.Temperature = atoi(stringvars[i].c_str());
    i++;
    //enable_recalc
    if(stringvars[i].compare("true")==0){
        params.Enable_recalc = true;
    }
    else if(stringvars[i].compare("false")==0){
        params.Enable_recalc = false;
    }
    else{
        cout << "Error setting event recalculation options" << endl;
        return false;
    }
    i++;
    params.Recalc_cutoff = atoi(stringvars[i].c_str());
    i++;
    //Tests
    //enable_exciton_diffusion_test
    if(stringvars[i].compare("true")==0){
        params.Enable_diffusion_test = true;
    }
    else if(stringvars[i].compare("false")==0){
        params.Enable_diffusion_test = false;
    }
    else{
        cout << "Error setting exciton diffusion test options" << endl;
        return false;
    }
    i++;
    params.N_tests = atoi(stringvars[i].c_str());
    i++;
    // Exciton Parameters
    params.Exciton_generation_rate = atof(stringvars[i].c_str());
    i++;
    params.Exciton_lifetime = atof(stringvars[i].c_str());
    i++;
    params.R_exciton_hopping = atof(stringvars[i].c_str());
    i++;
    params.FRET_cutoff = atoi(stringvars[i].c_str());
    i++;
    // Energetic Disorder Parameters
    //enable_gaussian_dos
    if(stringvars[i].compare("true")==0){
        params.Enable_gaussian_dos = true;
    }
    else if(stringvars[i].compare("false")==0){
        params.Enable_gaussian_dos = false;
    }
    else{
        cout << "Error setting Gaussian DOS options" << endl;
        return false;
    }
    i++;
    params.Site_energy_stdev = atof(stringvars[i].c_str());
    i++;
    //enable_exponential_dos
    if(stringvars[i].compare("true")==0){
        params.Enable_exponential_dos = true;
    }
    else if(stringvars[i].compare("false")==0){
        params.Enable_exponential_dos = false;
    }
    else{
        cout << "Error setting Exponential DOS options" << endl;
        return false;
    }
    i++;
    params.Site_energy_urbach = atof(stringvars[i].c_str());
    i++;
    // Error checking
    return true;
}


