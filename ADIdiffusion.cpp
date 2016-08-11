//
// Created by Isabelle Tan on 04-08-16.
//

#include "ADIdiffusion.h"

using namespace std;

// A function that solves a linear system in O(n), given the matrix is tridiagonal with a, b and c the diagonal constants
void ThomasAlg(const int N, const value_type a, const value_type b, const value_type c, Ref<MatrixXd> x, Ref<MatrixXd> r){
    // Allocate arrays for new constants
    value_type * y = new value_type[N-1];
    value_type * rho = new value_type[N];

    // FORWARD ELIMINATION
    // First row
    y[0] = c/b;
    rho[0] = r(0,0)/b;

    // Middle rows
    for (int i = 1; i < N-1; ++i) {
        y[i] = c/(b - a*y[i-1]);
        rho[i] = (r(i,0) - a*rho[i-1])/(b - a*y[i-1]);
    }

    // Last row
    rho[N-1] = (r(N-1,0) - a*rho[N-2])/(b - a*y[N-2]);

    // BACKWARD ELIMINATION
    // Last row
    x(N-1,0) = rho[N-1];

    // Rest of the rows
    for (int j = N-2; j >= 0; --j) {
        x(j,0) = rho[j] - y[j]*x(j+1,0);
    }
    
    // Free memory
    delete[] y;
    delete[] rho;

    return;
}

// This function computes one time iteration of diffusion using the ADI method.
// q_0 and q_new contain the data transposed
void ADI(Ref<MatrixXd> q_0, Ref<MatrixXd> q_new, const value_type dt, const value_type dx, const value_type v){
    // Preparations
    const value_type r = v * dt/(2*dx*dx); // r is the constant in the ADI discretization
    const int n = q_0.rows();

    // FIRST STAGE
    // Compute right hand side
    MatrixXd RHS(n,n);
    RHS.setZero(n,n);


    // Create matrix with (1-2r) on diagonal and r on sub- and superdiagonal
    MatrixXd rightMat(n,n);
    rightMat.setZero(n,n);
    rightMat(0,0) = 1-2*r;
    for (int i = 1; i < n; ++i) {
        rightMat(i,i) = 1-2*r;
        rightMat(i, i-1) = r;
        rightMat(i-1,i) = r;
    }


    // multiply to obtain the RHS
    RHS = rightMat*q_0.transpose();


    // Solve left hand side with Thomas algorithm (per row)
    // Make variables containing the subdiagonal, diagonal and superdiagonal values
    const value_type a = -r;
    const value_type b = 1+2*r;
    const value_type c = -r;
    MatrixXd q_half(n,n);
    q_half.setZero(n,n);

    // Loop over columns in q_half
    for (int i = 0; i < n; ++i) {
        ThomasAlg(n, a, b, c, q_half.col(i), RHS.col(i));
    }


    // SECOND STAGE
    // Compute right hand side
    RHS = rightMat*q_half.transpose();

    // Solve left hand side with Thomas algorithm to compute the tranposed of the result (can't compute result directly
    // because Eigen can only take .col() and not .row()
    MatrixXd q_new_T(n,n);
    for (int j = 0; j < n; ++j) {
        ThomasAlg(n, a, b, c, q_new_T.col(j), RHS.col(j));
    }


    // Take the transpose of the new data so it's columnwise again
    //q_new = q_new_T.transpose();
    q_new = q_new_T;
    return;
}