//
// Created by Isabelle Tan on 04-08-16.
//

#ifndef VORTEX_ADIDIFFUSION_TEST_H
#define VORTEX_ADIDIFFUSION_TEST_H

// Set the type to double
typedef double value_type;

// A function that solves a linear system in O(n), given the matrix is tridiagonal with a, b and c the diagonal constants
bool ThomasAlg_test();


// A function to reorder the particle vectors from row to column wise or the other way around (improves computation time
// in 2nd step of ADI diffusion method.)
bool changeOrdering_test();


// A function to perform one ADI diffusion time iteration
bool ADI_test(value_type * const q);


// A function to set the boundary values for the diffusion
bool setBoundaries_test(value_type * const q);


#endif //VORTEX_ADIDIFFUSION_TEST_H
