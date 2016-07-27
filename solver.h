//
// Created by Isabelle Tan on 27-07-16.
//

#ifndef VORTEX_SOLVER_H
#define VORTEX_SOLVER_H

// Set the value_type
typedef double value_type;

/*
 * This function computes the potential from all source particles at the target locations
 */
void potential();

/*
 * This function computes the velocity of the particles by taking the curl of the streamfunction.
 */
void velocity(const int N, const value_type h, value_type * const u, value_type * const v, value_type * const phi);

/*
 * This function computes the vorticity of the particles by taking the curl of the velocity.
 */
void vorticity();

/*
 * This function computes the spread of vorticity resulting from diffusion, using a Crank Nicholson scheme.
 */
void diffusion();

/*
 * This function computes the new particle locations from the velocity by performing one explicit Euler step.
 */
void advection();

/*
 * This function computes the time-step size for the next iteration based on the current velocities.
 */
void timeStep();

#endif //VORTEX_SOLVER_H
