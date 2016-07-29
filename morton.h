
#ifndef VORTEX_MORTON_H
#define VORTEX_MORTON_H

#include <algorithm>		
#include <parallel/algorithm>
#include <omp.h>
#include <stdlib.h>
#include <thread>

#include "simulation.h"		// because of the parameter "depth"

/*
 * Compute the extent, i.e. the size of smallest square including all data points
 */
void extent(const int N, const double* const x, const double* const y, double& xmin, double& ymin, double& ext); 								// set 2, question 1, a)

/*
 * Assign each particle its Morton-index
 */	
void morton(const int N, const double* const x, const double* const y, const double xmin, const double ymin, const double ext, int* index); 	// set 2, question 1, b)

/*
 * Sort the particles according t their Morton-index
 */
void sort(const int N, int* index, int* keys) ;																									// set 2, question 1, c)

/*
 * Reorder them
 */
void reorder(const int N, const int* const keys, const double* const x, const double* const y, double* xsorted, double* ysorted);	 			// set 2, question 1, d)

#endif //VORTEX_MORTON_H
