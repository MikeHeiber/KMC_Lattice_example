// Copyright (c) 2017-2019 Michael C. Heiber
// This source file is part of the KMC_Lattice_example project, which is subject to the MIT License.
// For more information, see the LICENSE file that accompanies this software.
// The KMC_Lattice_example project can be found on Github at https://github.com/MikeHeiber/KMC_Lattice_example

#include "Exciton_sim.h"
#include "Parameters.h"
#include "Utils.h"
#include <mpi.h>
#include <fstream>
#include <string>
#include <vector>

using namespace std;
using namespace KMC_Lattice;
using namespace KMC_Lattice_example;

int main(int argc, char *argv[]) {
	// Initialize mpi options
	cout << "Initializing MPI options... ";
	MPI_Init(&argc, &argv);
	int nproc = 1;
	int procid = 0;
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	MPI_Comm_rank(MPI_COMM_WORLD, &procid);
	cout << "MPI initialization complete!" << endl;
	// Start timer
	auto time_start = time(NULL);
	// Import parameters and options from file and command line arguments
	cout << "Loading input parameters from file... " << endl;
	ifstream parameterfile(argv[1], ifstream::in);
	if (!parameterfile) {
		cout << "Error loading parameter file.  Program will now exit." << endl;
		return 0;
	}
	Parameters params;
	bool success = params.importParameters(parameterfile);
	parameterfile.close();
	if (!success) {
		cout << "Error importing parameters from parameter file.  Program will now exit." << endl;
		return 0;
	}
	cout << "Parameter loading complete!" << endl;
	// Initialize Simulation
	cout << procid << ": Initializing simulation " << procid << "..." << endl;
	Exciton_sim sim(params, procid);
	cout << procid << ": Simulation initialization complete!" << endl;
	// Begin Simulation loop
	cout << procid << ": Starting simulation..." << endl;
	bool End_sim = false;
	while (!End_sim) {
		success = sim.executeNextEvent();
		if (!success) {
			cout << procid << ": Event execution failed, simulation will now terminate." << endl;
			break;
		}
		// Check if simulation has finished
		End_sim = sim.checkFinished();
		// Output status
		if (sim.getN_events_executed() % 100000 == 0) {
			sim.outputStatus();
		}
	}
	cout << procid << ": Simulation finished." << endl;
	auto time_end = time(NULL);
	auto elapsedtime = difftime(time_end, time_start);
	// Output simulation results for each processor
	ofstream resultsfile("results" + to_string(procid) + ".txt");
	resultsfile << "KMC_Lattice_example Results:\n";
	resultsfile << "Calculation time elapsed is " << (double)elapsedtime / 60.0 << " minutes.\n";
	resultsfile << sim.getTime() << " seconds have been simulated.\n";
	resultsfile << sim.getN_events_executed() << " events have been executed.\n";
	resultsfile << sim.getN_excitons_created() << " excitons have been created.\n";
	if (params.Enable_diffusion_test) {
		resultsfile << "Exciton diffusion test results:\n";
		resultsfile << "Exciton diffusion length is " << sim.calculateDiffusionLength_avg() << " ± " << sim.calculateDiffusionLength_stdev() << " nm\n";
	}
	resultsfile << endl;
	resultsfile.close();
	// Output overall analysis results from all processors
	vector<double> diffusion_data;
	if (params.Enable_diffusion_test) {
		diffusion_data = MPI_gatherVectors(sim.getDiffusionData());
	}
	if (procid == 0) {
		ofstream analysisfile("analysis_summary.txt");
		analysisfile << "KMC_Lattice_example Results Summary:" << endl;
		analysisfile << nproc * sim.getN_excitons_recombined() << " total excitons tested." << endl;
		if (params.Enable_diffusion_test) {
			analysisfile << "Overall exciton diffusion test results:\n";
			analysisfile << "Exciton diffusion length is " << vector_avg(diffusion_data) << " ± " << vector_stdev(diffusion_data) << " nm\n";
		}
		analysisfile << endl;
		analysisfile.close();
	}
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	return 0;
}
