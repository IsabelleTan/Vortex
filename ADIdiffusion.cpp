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