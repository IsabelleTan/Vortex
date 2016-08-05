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


// A function to perform a diffusion loop and write the intermetidate results to a file
void ADI_test_output();


#endif //VORTEX_ADIDIFFUSION_TEST_H
