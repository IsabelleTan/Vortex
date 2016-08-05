//
// Created by Isabelle Tan on 04-08-16.
//

#include "ADIdiffusion.h"
#include <iostream>
#include <cmath>

using namespace std;

// A function that solves a linear system in O(n), given the matrix is tridiagonal with a, b and c the diagonal constants
void ThomasAlg(const int N, const value_type a, const value_type b, const value_type c, value_type * const x, value_type * const r){
    // Allocate arrays for new constants
    value_type * y = new value_type[N-1];
    value_type * rho = new value_type[N];

    // FORWARD ELIMINATION
    // First row
    y[0] = c/b;
    rho[0] = r[0]/b;

    // Middle rows
    for (int i = 1; i < N-1; ++i) {
        y[i] = c/(b - a*y[i-1]);
        rho[i] = (r[i] - a*rho[i-1])/(b - a*y[i-1]);
    }

    // Last row
    rho[N-1] = (r[N-1] - a*rho[N-2])/(b - a*y[N-2]);

    // BACKWARD ELIMINATION
    // Last row
    x[N-1] = rho[N-1];

    // Rest of the rows
    for (int j = N-2; j >= 0; --j) {
        x[j] = rho[j] - y[j]*x[j+1];
    }
    
    // Free memory
    delete[] y;
    delete[] r;

    return;
}


// A function to reorder the particle vectors from row to column wise or the other way around (improves computation time
// in 2nd step of ADI diffusion method.)
void changeOrdering(const int N, value_type * const q_in, value_type * const q_out){
    const int n = sqrt(N);

    // Check whether N is a perfect square
    if (n*n != N){
        cout << "N is not a perfect square. Cannot reorder." << endl;
        return;
    }
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            q_out[n*i + j] = q_in[n*j + i];
        }
    }

    return;
}

// This function computes one time iteration of diffusion using the ADI method.
void ADI(int N, value_type * const q_0, value_type * const q_new, const value_type dt, const value_type dx, const value_type v){
    // Preparations
    const value_type r = v * dt/(2*dx*dx); // r is the constant in the ADI discretization
    const int n = sqrt(N);

    // Check whether N is a perfect square
    if(n*n != N){
        cout << "N is not a perfect square so the grid is not square. Not performing diffusion." << endl;
        return;
    }

    // FIRST STAGE
    // Compute right hand side
    value_type * RHS = new value_type[N];

    // Loop over the grid locations and compute finite differences for Y direction
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            // Check if this is the top boundary
            if(i == 0){
                RHS[i * n + j] = (1 - 2 * r) * q_0[i * n + j] + r*q_0[(i + 1) * n + j];
            } // Check if this is the bottom boundary
            else if (i==n-1){
                RHS[i * n + j] = (1 - 2 * r) * q_0[i * n + j] + r*q_0[(i - 1) * n + j];
            } // Interior values
            else {
                // Compute the right hand side
                RHS[i * n + j] = (1 - 2 * r) * q_0[i * n + j] + r * (q_0[(i + 1) * n + j] + q_0[(i - 1) * n + j]);
            }
        }
    }


    // Solve left hand side with Thomas algorithm
    // Make variables containing the subdiagonal, diagonal and superdiagonal values
    const value_type a = -r;
    const value_type b = 1+2*r;
    const value_type c = -r;
    value_type * q_half_r = new value_type[N];

    ThomasAlg(N, a, b, c, q_half_r, RHS);


    // SECOND STAGE
    // Reorder the vorticity vector such that it is ordered columnwise
    value_type * q_half_c = new value_type[N];
    changeOrdering(N, q_half_r, q_half_c);


    // Compute right hand side
    // Loop over the grid locations and compute finite differences in the X direction
    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < n; ++i) {
            // Check if this is the left boundary
            if(j == 0 ){
                RHS[j*n + i] = (1-2*r) * q_half_c[j*n+i] + r*q_half_c[(j+1)*n + i];
            } // Check if this is the right boundary
            else if(j==n-1){
                RHS[j*n + i] = (1-2*r) * q_half_c[j*n+i] + r*q_half_c[(j-1)*n + i];
            } // It must be an interior value
            else{
                // Compute the right hand side (re-use RHS vector, now ordered column wise)
                RHS[j*n + i] = (1-2*r) * q_half_c[j*n+i] + r * (q_half_c[(j-1)*n + i] + q_half_c[(j+1)*n + i]);
            }
        }
    }

    // Solve left hand side with Thomas algorithm
    value_type * q_new_c = new value_type[N];
    ThomasAlg(N, a, b, c, q_new_c, RHS);

    // Reorder the vector again so it is ordered row wise again
    changeOrdering(N, q_new_c, q_new);

    // Free allocated memory
    delete[] q_half_r;
    delete[] q_new_c;
    delete[] RHS;
    delete[] q_half_c;

    return;
}