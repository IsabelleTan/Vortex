//
// Created by Isabelle Tan on 04-08-16.
//

#ifndef VORTEX_ADIDIFFUSION_H
#define VORTEX_ADIDIFFUSION_H

#include <Eigen/Dense>

using namespace Eigen;

// Set the type to double
typedef double value_type;

// A function that solves a linear system in O(n), given the matrix is tridiagonal with a, b and c vectors containing the diagonal constants
void ThomasAlg(const int N, const value_type a, const value_type b, const value_type c, Ref<MatrixXd> x, Ref<MatrixXd> r);


// A function to perform one ADI diffusion time iteration
void ADI(Ref<MatrixXd> q_0, Ref<MatrixXd> q_new, const value_type dt, const value_type dx, const value_type v);


#endif //VORTEX_ADIDIFFUSION_H
