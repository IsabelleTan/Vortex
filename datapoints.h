#include <iostream> 
#include <iomanip>  	
#include <string>
#include <bitset>
#include <cstdio>

void fill_fromfile(const char* fname, int& N, float** x, float** y)
{
	// INPUT
	// fname 	name of file from which we want to fill vectors 
	// N		"empty"
	// x		non-memory-aligned empty vector of x-coordinates 
	// y		.................................. y-coordinates 
	// OUTPUT
	// N		# of data points (determined by file, not given to function as input) 
	// x		memory-aligned array of x-coordinates of data points read from file
	// y		....................... y-coordinates .............................
	
	FILE* f = fopen(fname, "r");			// open file to read from
	
	// determine n 
	char dummy[100]; 
	fscanf(f, "%s", dummy);					// read up to first whitespace? 
	fscanf(f, "%d", &N);					// determine the number of data points
	
	// do the memory-alignment 
	posix_memalign( (void**) x, 32, sizeof(float)*N );
	posix_memalign( (void**) y, 32, sizeof(float)*N );
	
	// fill up the vectors 
	fread(*x, sizeof(float), N, f);
	fread(*y, sizeof(float), N, f);
	
	// close the file stream
	fclose(f); 
}

void fill_random(int N, float** x, float** y)
{
	// INPUT 
	// N		# of data points
	// x		non memory-aligned empty array of x-coordinates 
	// y		................................. y-coordinates 
	// OUTPUT 
	// x		memory-aligned array filled with N random x-coordinates
	// y		......................................... y-coordinates
	
	// do the memory-alignment 
	posix_memalign( (void**) x, 32, sizeof(float)*N );
	posix_memalign( (void**) y, 32, sizeof(float)*N );
	
	// fill up the vectors 
	for(unsigned int i(0); i<N; i++){
		(*x)[i] = (drand48() - 0.5) * 10.; 	// drand48() returns a number between [0.0, 1.0)
		(*y)[i] = (drand48() - 0.5) * 10.; 	// so shift and expand the interval generated to [-5.0, 5.0).
	}
}


