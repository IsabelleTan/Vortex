//
// Created by Isabelle Tan on 27-07-16.
//

#ifndef VORTEX_FIELDS_H
#define VORTEX_FIELDS_H

#include <Eigen/Dense>

// Set the value_type
typedef double value_type;
using namespace Eigen;

/*
 * This function creates a 2D grid of x and y coordinates with domain (xRange,yRange) centered around 0
 */
bool grid(const int N, value_type * const x, value_type * const y, const value_type xRange, const value_type yRange);

/*
 * This function creates a Lamb-Oseen vortex.
 */
void lambOseen(const int N, value_type * const x, value_type * const y, value_type * const q, const value_type visc, const value_type circ, value_type t, const value_type xCenter, const value_type yCenter);

/*
 * A function to create a grid with two Lamb Oseen vortices. Both have coreRadius sigma=1 and are centered at (x1,y1) and (x2,y2).
 */
void lambOseenDouble(const int N, value_type * const x, value_type * const y, value_type * const q, const value_type visc, const value_type circ1, const value_type circ2, value_type t, const value_type x1, const value_type y1, const value_type x2, const value_type y2);

/*
 * A function to smoothly cut-off some value array after a certain radius, using a cosine function
 */
void smoothCutOff(const int N, value_type * const q_in, value_type * const x_in, value_type * const y_in, const value_type radius_start, const value_type radius_end);


/*
 * This function copies the content of an array x into an Eigen MatrixXd.
 * Ordering:
 * Array - lexicographical
 * Matrix - filled column wise (so transpose of what data locations look like)
 */
void arrayToMatrix(value_type * x, MatrixXd & M, bool invertOrder);

/*
 * This function copies the content of an Eigen MatrixXd into an array.
 * Ordering:
 * Array - lexicographical
 * Matrix - filled column wise (so transpose of what data locations look like)
 */
void matrixToArray(value_type * x, MatrixXd & M, bool invertOrder);

/*
 * This function computes the analytical solution and write the result to files
 */
void analyticalSolution();

#endif //VORTEX_FIELDS_H
