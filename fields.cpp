//
// Created by Isabelle Tan on 27-07-16.
//

#include "fields.h"
#include <iostream>



/*
 * This function creates a regular 2D grid of N (x,y) coordinates with domain (xRange,yRange) centered around 0. If it
 * succeeds it returns true, if it fails because N is not a perfect square it returns false.
 */
bool grid(const int N, value_type * const x, value_type * const y, const value_type xRange, const value_type yRange){
    // Compute the number of gridpoints in one row or column
    int M = (int)sqrt(N);

    // Test if M is an integer
    if (N != M*M){
        cout << "The number of gridpoints is not a perfect square, use a different N." << endl;
        return false;
    }

    // Compute parameters
    const value_type dx = xRange/(M-1);
    const value_type dy = yRange/(M-1);
    const value_type xLim = xRange/2;
    const value_type yLim = yRange/2;

    // Loop over y coordinates
    for (int i = 0; i < M; ++i) {
        // Loop over x coordinates
        for (int j = 0; j < M; ++j) {
            x[i*M + j] = j*dx - xLim;
            y[i*M + j] = -i*dy + yLim;

        }
    }

    return true;
}


/*
 * This function creates a Lamb-Oseen vortex.
 */
void lambOseen(const int N, value_type * const x, value_type * const y, value_type * const q, const value_type coreRadius, const value_type circ){
    // Compute the number of grid points in one row or column
    const int M = sqrt(N);

    // Assuming lexicographical ordering and the domain centered around (0,0)
    for (int i = 0; i < N; ++i) {
        // Compute the radius of the particle from the center of the domain
        value_type r = sqrt(x[i]*x[i] + y[i]*y[i]);
        q[i] = circ/(M_PI * coreRadius * coreRadius) * exp(-(r*r)/(coreRadius*coreRadius));
    }
}


/*
 * This function copies the content of an array x into an Eigen MatrixXd.
 * Ordering:
 * Array - lexicographical
 * Matrix - filled column wise (so transpose of what data locations look like)
 */
void arrayToMatrix(value_type * x, MatrixXd & M){
    // Loop over rows and columns
    const int r = M.rows();
    const int c = M.cols();
    for (int i = 0; i < r; ++i) {
        for (int j = 0; j < c; ++j) {
            M(j,i) = x[i*c + j];
        }
    }
    return;
}


/*
 * This function copies the content of an Eigen MatrixXd into an array.
 * Ordering:
 * Array - lexicographical
 * Matrix - filled column wise (so transpose of what data locations look like)
 */
void matrixToArray(value_type * x, MatrixXd & M){
    // Loop over rows and columns
    const int r = M.rows();
    const int c = M.cols();
    for (int i = 0; i < r; ++i) {
        for (int j = 0; j < c; ++j) {
            x[i*c + j] = M(j,i);
        }
    }
    return;
}