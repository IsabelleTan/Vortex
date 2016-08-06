//
// Created by Isabelle Tan on 27-07-16.
//

#ifndef VORTEX_FIELDS_H
#define VORTEX_FIELDS_H

#include <Eigen/Dense>

// Set the value_type
typedef double value_type;
using namespace Eigen;
using namespace std;

/*
 * This function creates a 2D grid of x and y coordinates with domain (xRange,yRange) centered around 0
 */
bool grid(const int N, value_type * const x, value_type * const y, const value_type xRange, const value_type yRange);

/*
 * This function creates a Lamb-Oseen vortex.
 */
void lambOseen(const int N, value_type * const x, value_type * const y, value_type * const q, const value_type coreRadius, const value_type circ);

/*
 * This function copies the content of an array x into an Eigen MatrixXd.
 * Ordering:
 * Array - lexicographical
 * Matrix - filled column wise (so transpose of what data locations look like)
 */
void arrayToMatrix(value_type * x, MatrixXd & M);

/*
 * This function copies the content of an Eigen MatrixXd into an array.
 * Ordering:
 * Array - lexicographical
 * Matrix - filled column wise (so transpose of what data locations look like)
 */
void matrixToArray(value_type * x, MatrixXd & M);

#endif //VORTEX_FIELDS_H
