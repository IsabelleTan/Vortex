//
// Created by Isabelle Tan on 27-07-16.
//

#ifndef VORTEX_EXPANSION_H
#define VORTEX_EXPANSION_H

#include "simulation.h"			// for simulation-parameters etc.

/*
 * The particle to expansion kernel. This function computes the alpha_k coefficients and stores their real and imaginary
 * parts in the arrays rexpansion and iexpansion.
 */
void p2e(const value_type* const x, const value_type* const y, const value_type* const q, const int nsources, const int order, const value_type xcom, const value_type ycom, value_type* const rexpansion, value_type* const iexpansion);

/*
 * The expansion to target particle kernel. This function computes the potential (streamfunction) value at the target particle located
 * at xtarget, ytarget.
 */
value_type e2p(const double xtarget, const double ytarget, const double q, const int order, const double* const rxps, const double* const ixps);

/*
 * The particle to particle kernel. This function computes the interactions between nsources source particles and one target particle.
 */
value_type p2p(const value_type* const xsources, const value_type* const ysources, const value_type* const q, const int nsources, const value_type xtarget, const value_type ytarget);


#endif //VORTEX_EXPANSION_H
