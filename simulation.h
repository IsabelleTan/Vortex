
#ifndef SIMULATION_H
#define SIMULATION_H

/*
 * Define the parameters that will be needed throughout the simulation
 */

// PI
#define MPI 3.141593

// Particles
#define nParticles 40000
#define deltaX 0.1
#define viscosity 0.1

// Time
#define t_0 0
#define deltaT 0.0001
#define timeIterations 25
#define writeFreq 1

// Initial Condition
#define coreRadius 1
#define circulation 10
#define xCenter1 0
#define yCenter1 0

// Cut Off
#define startRatio 0.6
#define endRatio 0.9

// Quadtree
#define depthtree 16 		// depth of quadtree // choice of levels: Beatson Greengard suggest approx log_2(N), where N = number of particles
#define kleaf 32			// leaf capacity

// Multipole
#define exp_order 12		// expansion order of the multipole (determined by convergence analysis)
#define theta_dist 0.5		// used to control which nodes are considered as "far"

typedef double value_type;

/*
 * Run the simulation 
 */
void run_simulation();

/*
 * Time the simulation
 */
void time_simulation(int nPart, int nSim);



#endif // SIMULATION_H
