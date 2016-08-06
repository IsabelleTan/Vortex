//
// Created by Isabelle Tan on 04-08-16.
//

#ifndef VORTEX_ADIDIFFUSION_TEST_H
#define VORTEX_ADIDIFFUSION_TEST_H



// A function that solves a linear system in O(n), given the matrix is tridiagonal with a, b and c the diagonal constants
bool ThomasAlg_test();


// A function to perform a diffusion loop and write the intermetidate results to a file
void ADI_test_output();


#endif //VORTEX_ADIDIFFUSION_TEST_H
