
#ifndef VORTEX_MORTON_H
#define VORTEX_MORTON_H

#include <algorithm>		
#include <parallel/algorithm>
#include <omp.h>
#include <stdlib.h>
#include <thread>

#include "simulation.h"		// because of the parameter "depth"

/*
 * Put some random values in the arrays
 */
void initialize(const int N, value_type* x, value_type* y, value_type* mass);

/*
 * Compute the extent, i.e. the size of smallest square including all data points
 */
void extent(const int N, const double* const x, const double* const y, double& xmin, double& ymin, double& ext); 								// set 2, question 1, a)

/*
 * Assign each particle its Morton-index
 */	
void morton(const int N, const double* const x, const double* const y, const double xmin, const double ymin, const double ext, unsigned int* index); 	// set 2, question 1, b)

/*
 * Sort the particles according t their Morton-index
 */
void sort(const int N, unsigned int* index, int* keys) ;																									// set 2, question 1, c)

/*
 * Reorder them
 */
void reorder(const int N, const int* const keys, const double* const x, const double* const y, const double* const q, double* xsorted, double* ysorted, double* qsorted);// set 2, question 1, d)

#endif //VORTEX_MORTON_H
