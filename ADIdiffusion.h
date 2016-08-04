//
// Created by Isabelle Tan on 04-08-16.
//

#ifndef VORTEX_ADIDIFFUSION_H
#define VORTEX_ADIDIFFUSION_H

// Set the type to double
typedef double value_type;

// A function that solves a linear system in O(n), given the matrix is tridiagonal with a, b and c the diagonal constants
void ThomasAlg(const int N, const value_type a, const value_type b, const value_type c, value_type * const x, value_type * const r);


// A function to reorder the particle vectors from row to column wise or the other way around (improves computation time
// in 2nd step of ADI diffusion method.)
void changeOrdering(value_type * const q_in, value_type * const q_out);


// A function to perform one ADI diffusion time iteration
void ADI(value_type * const q);


// A function to set the boundary values for the diffusion
void setBoundaries(value_type * const q);

#endif //VORTEX_ADIDIFFUSION_H
