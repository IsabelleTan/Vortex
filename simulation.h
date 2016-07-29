
#ifndef SIMULATION_H
#define SIMULATION_H

/*
 * Define the parameters that will be needed throughout the simulation
 */

#define theta 0.5			// used to control which nodes are considered as "far"
#define exp_order 10		// expansion order of the multipole
#define depth 15 			// depth of quadtree // choice of levels: Beatson Greengard suggest approx log_2(N), where N = number of particles
#define k 32				// leaf capacity

/*
 * Run the simulation 
 */
void run_simulation(); 

#endif // SIMULATION_H
