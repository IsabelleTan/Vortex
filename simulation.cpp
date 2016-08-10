
#include <iostream>
#include "simulation.h"
#include "fields.h"
#include "datapoints.h"
#include "multipole.h"
#include "solver.h"
#include "ADIdiffusion.h"

void run_simulation(){
	// INITIALIZATION
	// Allocate and align memory
	value_type * x_source;
	value_type * y_source;
	value_type * q_source;
	value_type * x_target;
	value_type * y_target;
	value_type * u_target;
	value_type * v_target;
	value_type * q_target;
	MatrixXd q_targetM(static_cast<int>(sqrt(nParticles)),static_cast<int>(sqrt(nParticles)));
	value_type * q_diffused;
	MatrixXd q_diffusedM;
	value_type * pot_target;
	posix_memalign((void **)&x_source, 32, sizeof(double) * nParticles);
	posix_memalign((void **)&y_source, 32, sizeof(double) * nParticles);
	posix_memalign((void **)&q_source, 32, sizeof(double) * nParticles);
	posix_memalign((void **)&x_target, 32, sizeof(value_type) * nParticles);
	posix_memalign((void **)&y_target, 32, sizeof(value_type) * nParticles);
	posix_memalign((void **)&u_target, 32, sizeof(value_type) * nParticles);
	posix_memalign((void **)&v_target, 32, sizeof(value_type) * nParticles);
	posix_memalign((void **)&q_target, 32, sizeof(value_type) * nParticles);
	posix_memalign((void **)&q_diffused, 32, sizeof(value_type)*nParticles);
	posix_memalign((void **)&pot_target, 32, sizeof(value_type) * nParticles);

	// Initialize values with Lamb Oseen vortex
	grid(nParticles, x_source, y_source, deltaX*nParticles, deltaX*nParticles);
	lambOseen(nParticles,x_source,y_source,q_source,viscosity,circulation,0);
	std::string filenameX = "0_X.txt";
	std::string filenameY = "0_Y.txt";
	std::string filenameQ = "0_Q.txt";
	// write_to_file(filenameX, nParticles, x_source);
	// write_to_file(filenameY, nParticles, y_source);
	// write_to_file(filenameQ, nParticles, q_source);

	// TIME ITERATIONS
	for (int i = 1; i <= timeIterations; ++i) {
		// Create target grid
		grid(nParticles, x_target, y_target, deltaX*nParticles, deltaX*nParticles);

		// Compute the potential at the target locations
		potential(theta_dist,x_source,y_source,q_source,nParticles,x_target,y_target,nParticles, pot_target);

		// Compute the velocity as the curl of the streamfunction potential
		velocity(nParticles, deltaX, u_target, v_target, pot_target);

		// Compute the vorticity as the curl of the velocity
		vorticity(nParticles, deltaX, u_target, v_target, q_target);

		// Perform one diffusion iteration (convert to and from matrices for the ADI solver)
//			std::cout << "hello 1, here is q_target:  " << q_target << std::endl ; 
		arrayToMatrix(q_target, q_targetM, true);
//					std::cout << "and here is q_targetM:  " << q_targetM << std::endl ; 
			std::cout << "hello 2 " << std::endl ; 
		ADI(q_targetM, q_diffusedM, deltaT, deltaX, viscosity);	
			std::cout << "hello 3 " << std::endl ; 
		matrixToArray(q_diffused, q_diffusedM, true);
			std::cout << "hello 4 " << std::endl ; 

		// Perform one advection iteration
		advection(nParticles, deltaT, u_target, v_target, x_target, y_target);


		// Check if for this iteration the output should be written to a file
		// If so; write output to file
		if(i%writeFreq == 0){
			filenameX = std::to_string(i) + "_X.txt";
			filenameY = std::to_string(i) + "_Y.txt";
			filenameQ = std::to_string(i) + "_Q.txt";
			//write_to_file(filenameX.c_str(), nParticles, x_target);
			//write_to_file(filenameY.c_str(), nParticles, y_target);
			//write_to_file(filenameQ.c_str(), nParticles, q_diffused);
		}


		// Set values for next iteration
		x_source = x_target;
		y_source = y_target;
		q_source = q_diffused;
	}

	delete[] x_source;
	delete[] y_source;
	delete[] q_source;
	delete[] x_target;
	delete[] y_target;
	delete[] u_target;
	delete[] v_target;
	delete[] q_target;
	delete[] q_diffused;
	delete[] pot_target;

	return;
	
}

void time_simulation(int nPart, int nSim){
	// TODO finish this function that times every step of the simulation

	value_type totalTime = 0;
	value_type gridTime = 0;
	value_type potentialTime = 0;
	value_type velocityTime = 0;
	value_type vorticityTime = 0;
	value_type diffusionTime = 0;
	value_type advectionTime = 0;

	// Loop over number of simulations
	for (int i = 0; i < nSim; ++i) {
		// Run simulation and time every step

		// Add elapsed times to timing variables

	}

	// Compute averages and percentages
	totalTime/=nSim;
	gridTime/=nSim;
	potentialTime/=nSim;
	velocityTime/=nSim;
	vorticityTime/=nSim;
	diffusionTime/=nSim;
	advectionTime/=nSim;


	std::cout << "Finished performing " << nSim << " simulations. \n" << "Total time: " << totalTime << std::endl;
	std::cout << "Computation time kernels: " << std::endl;
	std::cout << "Grid: \t\t Time = " << gridTime << "\t\t % of total = " << gridTime/totalTime * 100 << std::endl;
	std::cout << "Potential: \t\t Time = " << potentialTime << "\t\t % of total = " << potentialTime/totalTime * 100 << std::endl;
	std::cout << "Velocity: \t\t Time = " << velocityTime << "\t\t % of total = " << velocityTime/totalTime * 100 << std::endl;
	std::cout << "Vorticity: \t\t Time = " << vorticityTime << "\t\t % of total = " << vorticityTime/totalTime * 100 << std::endl;
	std::cout << "Diffusion: \t\t Time = " << diffusionTime << "\t\t % of total = " << diffusionTime/totalTime * 100 << std::endl;
	std::cout << "Advection: \t\t Time = " << advectionTime << "\t\t % of total = " << advectionTime/totalTime * 100 << std::endl;

}

int main(){
	run_simulation();
	std::cout <<"Finished!" << std::endl;
	return 0;
}

/*
 * Overview of steps


1. Initialize particles with positions and masses.

2.


 */
