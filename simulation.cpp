
#include <iostream>
#include "simulation.h"
#include "fields.h"
#include "datapoints.h"
#include "multipole.h"
#include "solver.h"
#include "ADIdiffusion.h"

void run_simulation(){
	// INITIALIZATION

	// make sure the number of particles is an integer squared, to avoid problems in matrices
	int dim = static_cast<int>(sqrt(nParticles));
	assert(dim*dim == nParticles);

	value_type domainRange = (dim-1)*deltaX;

	// Allocate and align memory
	value_type * x_source;
	value_type * y_source;
	value_type * q_source;
	value_type * rhs;
	value_type * x_target;
	value_type * y_target;
	value_type * u_target;
	value_type * v_target;
	value_type * q_target;
	MatrixXd q_targetM(dim,dim);
	value_type * q_diffused;
	MatrixXd q_diffusedM(dim,dim);
	value_type * pot_target;
	posix_memalign((void **)&x_source, 32, sizeof(value_type) * nParticles);
	posix_memalign((void **)&y_source, 32, sizeof(value_type) * nParticles);
	posix_memalign((void **)&q_source, 32, sizeof(value_type) * nParticles);
	posix_memalign((void **)&rhs, 32, sizeof(value_type) * nParticles);
	posix_memalign((void **)&x_target, 32, sizeof(value_type) * nParticles);
	posix_memalign((void **)&y_target, 32, sizeof(value_type) * nParticles);
	posix_memalign((void **)&u_target, 32, sizeof(value_type) * nParticles);
	posix_memalign((void **)&v_target, 32, sizeof(value_type) * nParticles);
	posix_memalign((void **)&q_target, 32, sizeof(value_type) * nParticles);
	posix_memalign((void **)&q_diffused, 32, sizeof(value_type)*nParticles);
	posix_memalign((void **)&pot_target, 32, sizeof(value_type) * nParticles);

//! ---------------------------------------------------------- initial conditions = 0
	// Initialize values with Lamb Oseen vortex
	grid(nParticles, x_source, y_source, deltaX*(dim-1), deltaX*(dim-1));
	lambOseen(nParticles,x_source,y_source,q_source,viscosity,circulation,0);
	std::string filenameX = "0_X.txt";
	std::string filenameY = "0_Y.txt";
	std::string filenameQ = "0_Q.txt";
	std::string filenameV = "0_V.txt";
	std::string filenameP = "0_P.txt";
	write_to_file(filenameX.c_str(), nParticles, x_source);
	write_to_file(filenameY.c_str(), nParticles, y_source);
	write_to_file(filenameQ.c_str(), nParticles, q_source);
	write_to_file(filenameP.c_str(), nParticles, pot_target);
	value_type * vel;
	posix_memalign((void **)&vel, 32, sizeof(value_type) * nParticles);
	for(size_t i(0); i<nParticles; ++i){
		vel[i] = sqrt(v_target[i]*v_target[i] + u_target[i]*u_target[i]);
	}
	write_to_file(filenameV.c_str(), nParticles, vel);

	// TIME ITERATIONS
	for (int i = 1; i <= 1 /*timeIterations*/; ++i) {

		std::cout << "TIME iteration # " << i << std::endl;

		// Create target grid
		grid(nParticles, x_target, y_target, deltaX*(dim-1), deltaX*(dim-1));

//! ---------------------------------------------------------- after grid = 1 
		filenameX = std::to_string(1) + "_X.txt";
		filenameY = std::to_string(1) + "_Y.txt";
		filenameQ = std::to_string(1) + "_Q.txt";
		filenameV = std::to_string(1) + "_V.txt";
		filenameP = std::to_string(1) + "_P.txt";
		write_to_file(filenameP.c_str(), nParticles, pot_target);
		write_to_file(filenameX.c_str(), nParticles, x_target);
		write_to_file(filenameY.c_str(), nParticles, y_target);
		write_to_file(filenameQ.c_str(), nParticles, q_source);
// compute scalar velocity
		for(size_t i(0); i<nParticles; ++i){
			vel[i] = sqrt(v_target[i]*v_target[i] + u_target[i]*u_target[i]);
		}
		write_to_file(filenameV.c_str(), nParticles, vel);

		// Compute the potential at the target locations
		for(size_t i(0); i<nParticles; ++i){
			rhs[i] = -0.5 * (1./MPI) * deltaX * deltaX * q_source[i] ;
		}

		potential(theta_dist,x_source,y_source,rhs,nParticles,x_target,y_target,nParticles, pot_target);

//! ---------------------------------------------------------- after potential = 2		
		filenameX = std::to_string(2) + "_X.txt";
		filenameY = std::to_string(2) + "_Y.txt";
		filenameQ = std::to_string(2) + "_Q.txt";
		filenameV = std::to_string(2) + "_V.txt";
		filenameP = std::to_string(2) + "_P.txt";
		write_to_file(filenameP.c_str(), nParticles, pot_target);
		write_to_file(filenameX.c_str(), nParticles, x_target);
		write_to_file(filenameY.c_str(), nParticles, y_target);
		write_to_file(filenameQ.c_str(), nParticles, q_source);
// compute scalar velocity
		for(size_t i(0); i<nParticles; ++i){
			vel[i] = sqrt(v_target[i]*v_target[i] + u_target[i]*u_target[i]);
		}
		write_to_file(filenameV.c_str(), nParticles, vel);

		// Compute the velocity as the curl of the streamfunction potential
		velocity(nParticles, deltaX, u_target, v_target, pot_target);											////////////// ****

/*for(size_t i=0; i<nParticles; i++){
	std::cout << "index i = " << i << std::endl; 
	assert( (std::isfinite(u_target[i])) && (std::isfinite(v_target[i])) );
}*/


//! ---------------------------------------------------------- after velocity = 3 		
		filenameX = std::to_string(3) + "_X.txt";
		filenameY = std::to_string(3) + "_Y.txt";
		filenameQ = std::to_string(3) + "_Q.txt";
		filenameV = std::to_string(3) + "_V.txt";
		filenameP = std::to_string(3) + "_P.txt";
		write_to_file(filenameP.c_str(), nParticles, pot_target);
		write_to_file(filenameX.c_str(), nParticles, x_target);
		write_to_file(filenameY.c_str(), nParticles, y_target);
		write_to_file(filenameQ.c_str(), nParticles, q_source);
// compute scalar velocity
		for(size_t i(0); i<nParticles; ++i){
			vel[i] = sqrt(v_target[i]*v_target[i] + u_target[i]*u_target[i]);
		}
		write_to_file(filenameV.c_str(), nParticles, vel);

		// Compute the vorticity as the curl of the velocity
		vorticity(nParticles, deltaX, u_target, v_target, q_target);

//! ---------------------------------------------------------- after vorticity = 4

		filenameX = std::to_string(4) + "_X.txt";
		filenameY = std::to_string(4) + "_Y.txt";
		filenameQ = std::to_string(4) + "_Q.txt";
		filenameV = std::to_string(4) + "_V.txt";
		filenameP = std::to_string(4) + "_P.txt";
		write_to_file(filenameP.c_str(), nParticles, pot_target);
		write_to_file(filenameX.c_str(), nParticles, x_target);
		write_to_file(filenameY.c_str(), nParticles, y_target);
		write_to_file(filenameQ.c_str(), nParticles, q_target);

		// Perform a smooth cut off of the vorticity
		smoothCutOff(nParticles, q_target, x_target, y_target, startRatio*(domainRange/2), endRatio*(domainRange/2));

//! ---------------------------------------------------------- after cut off = 5

		filenameX = std::to_string(5) + "_X.txt";
		filenameY = std::to_string(5) + "_Y.txt";
		filenameQ = std::to_string(5) + "_Q.txt";
		filenameV = std::to_string(5) + "_V.txt";
		filenameP = std::to_string(5) + "_P.txt";
		write_to_file(filenameP.c_str(), nParticles, pot_target);
		write_to_file(filenameX.c_str(), nParticles, x_target);
		write_to_file(filenameY.c_str(), nParticles, y_target);
		write_to_file(filenameQ.c_str(), nParticles, q_target);
// compute scalar velocity
		for(size_t i(0); i<nParticles; ++i){
			vel[i] = sqrt(v_target[i]*v_target[i] + u_target[i]*u_target[i]);
		}
		write_to_file(filenameV.c_str(), nParticles, vel);

		// Perform one diffusion iteration (convert to and from matrices for the ADI solver)
		arrayToMatrix(q_target, q_targetM, true);
		ADI(q_targetM, q_diffusedM, deltaT, deltaX, viscosity);
		matrixToArray(q_diffused, q_diffusedM, true);

//! ---------------------------------------------------------- after ADI diff = 5		
		filenameX = std::to_string(6) + "_X.txt";
		filenameY = std::to_string(6) + "_Y.txt";
		filenameQ = std::to_string(6) + "_Q.txt";
		filenameV = std::to_string(6) + "_V.txt";
		filenameP = std::to_string(6) + "_P.txt";
		write_to_file(filenameP.c_str(), nParticles, pot_target);
		write_to_file(filenameX.c_str(), nParticles, x_target);
		write_to_file(filenameY.c_str(), nParticles, y_target);
		write_to_file(filenameQ.c_str(), nParticles, q_diffused);
// compute scalar velocity
		for(size_t i(0); i<nParticles; ++i){
			vel[i] = sqrt(v_target[i]*v_target[i] + u_target[i]*u_target[i]);
		}
		write_to_file(filenameV.c_str(), nParticles, vel);

		// Perform one advection iteration
		advection(nParticles, deltaT, u_target, v_target, x_target, y_target);

//! ---------------------------------------------------------- after advection = 6		
		filenameX = std::to_string(7) + "_X.txt";
		filenameY = std::to_string(7) + "_Y.txt";
		filenameQ = std::to_string(7) + "_Q.txt";
		filenameV = std::to_string(7) + "_V.txt";
		filenameP = std::to_string(7) + "_P.txt";
		write_to_file(filenameP.c_str(), nParticles, pot_target);
		write_to_file(filenameX.c_str(), nParticles, x_target);
		write_to_file(filenameY.c_str(), nParticles, y_target);
		write_to_file(filenameQ.c_str(), nParticles, q_diffused);
// compute scalar velocity
		for(size_t i(0); i<nParticles; ++i){
			vel[i] = sqrt(v_target[i]*v_target[i] + u_target[i]*u_target[i]);
		}
		write_to_file(filenameV.c_str(), nParticles, vel);
		free(vel);

/*		// Check if for this iteration the output should be written to a file
		// If so; write output to file
		if(i%writeFreq == 0){
			filenameX = std::to_string(i) + "_X.txt";
			filenameY = std::to_string(i) + "_Y.txt";
			filenameQ = std::to_string(i) + "_Q.txt";
			filenameV = std::to_string(i) + "_V.txt";
			write_to_file(filenameX.c_str(), nParticles, x_target);
			write_to_file(filenameY.c_str(), nParticles, y_target);
			write_to_file(filenameQ.c_str(), nParticles, q_diffused);
			// compute scalar velocity
			value_type * vel;
			posix_memalign((void **)&vel, 32, sizeof(value_type) * nParticles);
			for(size_t i(0); i<nParticles; ++i){
				vel[i] = sqrt(v_target[i]*v_target[i] + u_target[i]*u_target[i]);
			}
			write_to_file(filenameV.c_str(), nParticles, vel);
			free(vel);
		} */

		// Set values for next iteration
		x_source = x_target;
		y_source = y_target;
		q_source = q_diffused;

	}

	free(x_source);
	free(y_source);
	free(q_source);
//	free(x_target); should not be freed, since it has the same address as x_source which has already been freed
//	free(y_target); should not be freed, since it has the same address as y_source which has already been freed
	free(u_target);
	free(v_target);
	free(q_target);
//	free(q_diffused); should not be freed, since it has the same address as q_source which has already been freed
	free(pot_target);

	return;

}

// Use special size for nParticles!
void time_simulation(int nPart, int nSim){
	// TODO finish this function that times every step of the simulation
	// INITIALIZATION

	// Timing variables
	value_type totalTime = 0;
	value_type gridTime = 0;
	value_type potentialTime = 0;
	value_type velocityTime = 0;
	value_type vorticityTime = 0;
	value_type diffusionTime = 0;
	value_type advectionTime = 0;
	auto start = std::chrono::high_resolution_clock::now();
	auto end   = std::chrono::high_resolution_clock::now();

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
	posix_memalign((void **)&x_source, 32, sizeof(double) * nPart);
	posix_memalign((void **)&y_source, 32, sizeof(double) * nPart);
	posix_memalign((void **)&q_source, 32, sizeof(double) * nPart);
	posix_memalign((void **)&x_target, 32, sizeof(value_type) * nPart);
	posix_memalign((void **)&y_target, 32, sizeof(value_type) * nPart);
	posix_memalign((void **)&u_target, 32, sizeof(value_type) * nPart);
	posix_memalign((void **)&v_target, 32, sizeof(value_type) * nPart);
	posix_memalign((void **)&q_target, 32, sizeof(value_type) * nPart);
	posix_memalign((void **)&q_diffused, 32, sizeof(value_type) * nPart);
	posix_memalign((void **)&pot_target, 32, sizeof(value_type) * nPart);

	// Initialize values with Lamb Oseen vortex
	grid(nPart, x_source, y_source, deltaX*nPart, deltaX*nPart);
	lambOseen(nPart,x_source,y_source,q_source,viscosity,circulation,0);
	// Do not write to file during timing simulations
	//std::string filenameX = "0_X.txt";
	//std::string filenameY = "0_Y.txt";
	//std::string filenameQ = "0_Q.txt";
	// write_to_file(filenameX, nParticles, x_source);
	// write_to_file(filenameY, nParticles, y_source);
	// write_to_file(filenameQ, nParticles, q_source);
	// Add elapsed times to timing variables

	// TIME ITERATIONS
	for (int i = 1; i <= timeIterations; ++i) {
		// Create target grid
		start = std::chrono::high_resolution_clock::now();
		grid(nPart, x_target, y_target, deltaX*nPart, deltaX*nPart);
		end = std::chrono::high_resolution_clock::now();
		gridTime += static_cast<std::chrono::duration<double>>(end-start).count();
		totalTime += static_cast<std::chrono::duration<double>>(end-start).count();


		/*// Compute the potential at the target locations
		start = std::chrono::high_resolution_clock::now();
		potential(theta_dist,x_source,y_source,q_source,nPart,x_target,y_target,nPart, pot_target);
		end = std::chrono::high_resolution_clock::now();
		potentialTime += static_cast<std::chrono::duration<double>>(end-start).count();
		totalTime += static_cast<std::chrono::duration<double>>(end-start).count();

		// Compute the velocity as the curl of the streamfunction potential
		start = std::chrono::high_resolution_clock::now();
		velocity(nPart, deltaX, u_target, v_target, pot_target);
		end = std::chrono::high_resolution_clock::now();
		velocityTime += static_cast<std::chrono::duration<double>>(end-start).count();
		totalTime += static_cast<std::chrono::duration<double>>(end-start).count();

		// Compute the vorticity as the curl of the velocity
		start = std::chrono::high_resolution_clock::now();
		vorticity(nPart, deltaX, u_target, v_target, q_target);
		end = std::chrono::high_resolution_clock::now();
		vorticityTime += static_cast<std::chrono::duration<double>>(end-start).count();
		totalTime += static_cast<std::chrono::duration<double>>(end-start).count();

		// Perform one diffusion iteration (convert to and from matrices for the ADI solver)
		start = std::chrono::high_resolution_clock::now();
		arrayToMatrix(q_target, q_targetM, true);
		ADI(q_targetM, q_diffusedM, deltaT, deltaX, viscosity);
		matrixToArray(q_diffused, q_diffusedM, true);
		end = std::chrono::high_resolution_clock::now();
		diffusionTime += static_cast<std::chrono::duration<double>>(end-start).count();
		totalTime += static_cast<std::chrono::duration<double>>(end-start).count();

		// Perform one advection iteration
		start = std::chrono::high_resolution_clock::now();
		advection(nPart, deltaT, u_target, v_target, x_target, y_target);
		end = std::chrono::high_resolution_clock::now();
		advectionTime += static_cast<std::chrono::duration<double>>(end-start).count();
		totalTime += static_cast<std::chrono::duration<double>>(end-start).count(); */


		// Check if for this iteration the output should be written to a file
		// If so; write output to file
		/* Do not write to files during timing simulation runs
         * if(i%writeFreq == 0){
            filenameX = std::to_string(i) + "_X.txt";
            filenameY = std::to_string(i) + "_Y.txt";
            filenameQ = std::to_string(i) + "_Q.txt";
            write_to_file(filenameX.c_str(), nPart, x_target);
            write_to_file(filenameY.c_str(), nPart, y_target);
            write_to_file(filenameQ.c_str(), nPart, q_diffused);
        } */


		// Set values for next iteration
		x_source = x_target;
		y_source = y_target;
		q_source = q_diffused;
	}

	/* TODO If we want to simulate multiple times we could use put this loop around the code above. (Watch out with allocation and freeing!)
	// Loop over number of simulations
	for (int i = 0; i < nSim; ++i) {
		// Run simulation and time every step

	}

	// Compute averages and percentages
	totalTime/=nSim;
	gridTime/=nSim;
	potentialTime/=nSim;
	velocityTime/=nSim;
	vorticityTime/=nSim;
	diffusionTime/=nSim;
	advectionTime/=nSim;
	 */


	std::cout << "Finished performing " << nSim << " simulations with " << nPart << " particles. \n" << "Total time: " << totalTime << std::endl;
	std::cout << "Computation time kernels: " << std::endl;
	std::cout << "Grid: \t\t Time = " << gridTime << "\t\t % of total = " << gridTime/totalTime * 100 << std::endl;
	std::cout << "Potential: \t\t Time = " << potentialTime << "\t\t % of total = " << potentialTime/totalTime * 100 << std::endl;
	std::cout << "Velocity: \t\t Time = " << velocityTime << "\t\t % of total = " << velocityTime/totalTime * 100 << std::endl;
	std::cout << "Vorticity: \t\t Time = " << vorticityTime << "\t\t % of total = " << vorticityTime/totalTime * 100 << std::endl;
	std::cout << "Diffusion: \t\t Time = " << diffusionTime << "\t\t % of total = " << diffusionTime/totalTime * 100 << std::endl;
	std::cout << "Advection: \t\t Time = " << advectionTime << "\t\t % of total = " << advectionTime/totalTime * 100 << std::endl;


	// Free memory
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
	//time_simulation(100, 1);
	std::cout <<"Finished!" << std::endl;
	return 0;

}
