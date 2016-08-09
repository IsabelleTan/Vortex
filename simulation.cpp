
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
	MatrixXd q_targetM;
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
		arrayToMatrix(q_target, q_targetM, true);
		ADI(q_targetM, q_diffusedM, deltaT, deltaX, viscosity);
		matrixToArray(q_diffused, q_diffusedM, true);

		// Perform one advection iteration
		advection(nParticles, deltaT, u_target, v_target, x_target, y_target);


		// Check if for this iteration the output should be written to a file
		// If so; write output to file
		if(i%writeFreq == 0){
			filenameX = std::to_string(i) + "_X.txt";
			filenameY = std::to_string(i) + "_Y.txt";
			filenameQ = std::to_string(i) + "_Q.txt";
			write_to_file(filenameX.c_str(), nParticles, x_target);
			write_to_file(filenameY.c_str(), nParticles, y_target);
			write_to_file(filenameQ.c_str(), nParticles, q_diffused);
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

int main(){
	run_simulation();
	return 0;
}

/*
 * Overview of steps


1. Initialize particles with positions and masses.

2.


 */