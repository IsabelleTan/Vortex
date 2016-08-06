//
// Created by Isabelle Tan on 27-07-16.
//

#include "solver.h"
#include "cmath"

/*
 * This function computes the potential from all source particles at the target locations
 */
void potential(){
    // TODO complete potential
}

/*
 * This function computes the velocity of the particles by taking the curl of the streamfunction.
 * Central differences is used and for the boundary the zero gradient rule is used. The boundary corners are computed as
 * the average of their neighbours.
 */
void velocity(const int N, const value_type h, value_type * const u, value_type * const v, value_type * const phi){
    // Compute the number of particles in one row or column of the grid
    const int M = sqrt(N);

    // Assuming lexicographical ordering of the grid particles
    for (int i = 0; i < N; ++i) {
        // Check for boundary conditions
        if(i < M || i%M ==0 || i%M==M-1 || i > (M-1)*M){
            // Do nothing on the boundaries
        } else {
            // Compute the partial derivatives with central differences to obtain the u and v velocities
            u[i] = (phi[i-M] - phi[i+M])/(2*h);
            v[i] = (phi[i-1]-phi[i+1])/(2*h);
        }
    }

    // TODO I don't know about this
    // Set the boundary velocities equal to the inner ring velocities
    for (int j = 1; j < M-1; ++j) {
        // Top boundary
        u[j] = u[j+M];
        v[j] = v[j+M];
        // Left boundary
        u[j*M] = u[j*M + 1];
        v[j*M] = v[j*M +1];
        // Right boundary
        u[(j+1)*M - 1] = u[(j+1)*M -2];
        v[(j+1)*M - 1] = v[(j+1)*M -2];
        // Bottom boundary
        u[N-j-1] = u[N-j-1-M];
        v[N-j-1] = v[N-j-1-M];
    }

    // Set the boundary corners to the average of their neighbours
    u[0] = (u[1] + u[M])/2;
    u[M-1] = (u[M-2]+u[2*M-1])/2;
    u[N-M] = (u[N-M+1]+u[N-(2*M)])/2;
    u[N-1] = (u[N-2]+u[N-M-1])/2;

    v[0] = (v[1] + v[M])/2;
    v[M-1] = (v[M-2]+v[2*M-1])/2;
    v[N-M] = (v[N-M+1]+v[N-(2*M)])/2;
    v[N-1] = (v[N-2]+v[N-M-1])/2;
}

/*
 * This function computes the vorticity of the particles by taking the curl of the velocity.
 */
void vorticity(const int N, const value_type h, value_type * const u, value_type * const v, value_type * const q){
    // Compute the number of particles in one row or column
    const int n = sqrt(N);

    // Assuming lexicographical ordering of the grid particles
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            value_type u_y;
            value_type v_x;
            if (i == 0) {
                // Top boundary
                u_y = -u[(i + 1) * n + j] / (2 * h);
            } else if (i == n - 1) {
                // Bottom boundary
                u_y = u[(i - 1) * n + j] / (2 * h);
            } else {
                // Interior for y derivative
                u_y = (u[(i - 1) * n + j] - u[(i + 1) * n + j]) / (2 * h);
            }

            if (j == 0) {
                // Left boundary
                v_x = v[i * n + j + 1] / (2 * h);
            } else if (j == n - 1) {
                // Right boundary
                v_x = -v[i * n + j - 1] / (2 * h);
            } else {
                // Interior for x derivative
                v_x = (v[i * n + j + 1] - v[i * n + j - 1]) / (2 * h);
            }

            // Set the vorticity
            q[i*n+j] = v_x + u_y;
        }
    }
}

/*
 * This function computes the new particle locations from the velocity by performing one explicit Euler step.
 */
void advection(const int N, const value_type dt, value_type * const u, value_type * const v, value_type * const x, value_type * const y){
    // Use forward explicit Euler
    for (int i = 0; i < N; ++i) {
        x[i] += u[i] * dt;
        y[i] += v[i] * dt;
    }
}

/*
 * This function computes the time-step size for the next iteration based on the current velocities.
 */
void timeStep(double& dt){
    // TODO complete timeStep
    
    //! notes by Shoshana: 
    //! when we met the TAs at the end of the semester, they suggested we first simply use dt = 10^(-3) 
    //! as a time-step and once the simulation works, move on to an adaptive one.
    
    dt = 0.001; // seconds
}
