#ifndef VORTEX_DATAPOINTS_H
#define VORTEX_DATAPOINTS_H

#include <cstdio>
#include <stdlib.h>

typedef double value_type;

// load data from a binary file into arrays x and y 
// for layout, look at commentaries inside function body
void load_data_from_file(const char* fname, int& N, value_type** x, value_type** y);

// load N random floats into arrays x and y 
// the random values will be between [-5.0, 5.0)
void load_data_random(int N, value_type** x, value_type** y);


// clear a file of all its contents
// should be called in the beginning of a simulation to make sure we're working with clean files
// and aren't cluttered by data from previous simulations
void clean_file(const char* fname);

// write array x to file fname in binary
// if there is data already written in the file, new data is appended at the end of it
// for layout, look at commentaries inside function body
void write_to_file(const char* fname, int N, const value_type* x);

#endif // VORTEX_DATAPOINTS_H
