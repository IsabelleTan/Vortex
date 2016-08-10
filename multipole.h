
#ifndef VORTEX_MULTIPOLE_H
#define VORTEX_MULITPOLE_H

#include "simulation.h"
#include "expansion.h"
#include "quadtree.h"


/*
 * This function computes the potential from all source particles at the target locations.
 * This function computes the potential (i.e. the stream function psi) at all target locations (on a grid),
 * using the vorticity ("mass") of all source particles, which could be located anywhere in domain. 
 * We're using a recursive approach and starting with the root node
 */
//! parallel with respect to targets 
//! but load imbalance: traversal of tree is different for each target 
void potential(double theta, 
			  double *xsrc, double *ysrc, double *qsrc, int nsrc,
			  double *xdst, double *ydst, int ndst, double *potdst);


/*
 * Evaluate the potential (i.e. the stream function psi) at a particular target location (xt, yt), 
 * given an array of nodes and their expansions 
 */
void evaluate(const Node* nodes, const int node_id, const double *xdata, const double *ydata, const double *mdata,
		const double thetasquared, double * const result, const double xt, const double yt);

#endif //VORTEX_MULTIPOLE_H

