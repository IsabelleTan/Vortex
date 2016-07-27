//
// Created by Isabelle Tan on 27-07-16.
//

#ifndef VORTEX_FIELDS_H
#define VORTEX_FIELDS_H

// Set the value_type
typedef double value_type;

/*
 * This function creates a Lamb-Oseen vortex.
 */
void lambOseen(const int N, const value_type xMax, const value_type yMax, const value_type radius, const value_type coreRadius, const value_type circ, value_type * const q, const value_type nu);


#endif //VORTEX_FIELDS_H
