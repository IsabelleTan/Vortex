#include <iostream>
#include <iomanip>
#include <bitset>
#include "datapoints.h"
#include <unistd.h>
#define GetCurrentDir getcwd

// load data from a binary file into arrays x and y 
// for layout, look at commentaries inside function body
void load_data_from_file(const char* fname, int& N, value_type** x, value_type** y)
{
	// INPUT
	// fname 	name of file from which we want to read vectors points.
	//          note to user: give file name with its extension
	// N		"empty"
	// x		non-memory-aligned empty vector of x-coordinates 
	// y		.................................. y-coordinates 
	// OUTPUT
	// N		# of data points (determined by file, not given to function as input) 
	// x		memory-aligned array of x-coordinates of data points read from file
	// y		....................... y-coordinates .............................
	// DATA LAYOUT IN FILE
	// N (ENTER)						  (integer in decimal format) 
	// x[0]x[1]...x[N-1]y[0]y[1]...y[N-1] (where each point is a float written in binary)

	FILE* f = fopen(fname, "r");			// open file to read from in "read" mode

	// determine n
	char dummy[100];
	fscanf(f, "%s", dummy);					// read up to first whitespace? 
	fscanf(f, "%d\n", &N);					// read the number of data points in decimal format

	// do the memory-alignment 
	posix_memalign( (void**) x, 32, sizeof(value_type)*N );
	posix_memalign( (void**) y, 32, sizeof(value_type)*N );

	// fill up the vectors 
	fread(*x, sizeof(value_type), N, f);
	fread(*y, sizeof(value_type), N, f);

	// close the file stream
	fclose(f);
}

// load N random floats into arrays x and y 
// the random values will be between [-5.0, 5.0)
void load_data_random(int N, value_type** q)
{
	// INPUT 
	// N		# of data points
	// q		non memory-aligned empty array
	// OUTPUT 
	// q		memory-aligned array filled with N random values

	// do the memory-alignment 
	posix_memalign( (void**) q, 32, sizeof(value_type)*N );

	// fill up the vectors 
	for(unsigned int i(0); i<N; i++){
		(*q)[i] = (drand48() - 0.5) * 10.;
	}
}

// load N random floats into arrays x and y 
// the random values will be between [-5.0, 5.0)
void load_data_random(int N, value_type** x, value_type** y) {
	// INPUT 
	// N		# of data points
	// x		non memory-aligned empty array of x-coordinates 
	// y		................................. y-coordinates 
	// OUTPUT 
	// x		memory-aligned array filled with N random x-coordinates
	// y		......................................... y-coordinates

	// do the memory-alignment 
	posix_memalign((void **) x, 32, sizeof(value_type) * N);
	posix_memalign((void **) y, 32, sizeof(value_type) * N);

	// fill up the vectors 
	for (unsigned int i(0); i < N; i++) {
		(*x)[i] = (drand48() - 0.5) * 10.;    // drand48() returns a number between [0.0, 1.0)
		(*y)[i] = (drand48() - 0.5) * 10.;    // so shift and expand the interval generated to [-5.0, 5.0).
	}
}

// clear a file of all its contents
// should be called in the beginning of a simulation to make sure we're working with clean files
// and aren't cluttered by data from previous simulations
void clean_file(const char* fname){
	// INPUT
	// fname	name of file to be emptied
	//          note to user: give file name with its extension
	// OUTPUT
	// the file named "fname" is empty, but still exists after this function is called
	// to be used

	FILE* f = fopen(fname, "w"); 			// open the file in "write" mode
	// if a file with this name already exists, its contents are its contents are discarded and the file is treated as a new empty file.
	fclose(f);
}

// write array x to file fname in binary
// if there is data already written in the file, new data is appended at the end of it
// for layout, look at commentaries inside function body
void write_to_file(const char* fname, int N, const value_type* x){
    // INPUT
    // fname	name of file to be written in
    //          note to user: give file name with its extension
    // N		# of data points
    // x		non memory-aligned empty array of x-coordinates
    // DATA LAYOUT IN FILE
    // N (ENTER)						  (integer in decimal format)
    // x[0]x[1]...x[N-1]y[0]y[1]...y[N-1] (where each point is a float written in binary)
    // Note: if a file named "fname" already exists, this function appends data at the end of the file

    FILE* f = fopen(fname, "wb");			// open file in "write" mode
    if (f==NULL){
        std::cout << "Opening file " << fname << " failed." << std::endl;
    }

    // write N
    fprintf(f, "%d\n", N);					// write N in decimal format, then go to a new line ("\n")
    // write array content
    fwrite(x, sizeof(value_type), N, f);		// write the array contents in binary format

    std:: cout << "Written to file " << fname << std::endl;
    fclose(f);
}

void printWorkingDirectory(){
	char cCurrentPath[FILENAME_MAX];

	if (!GetCurrentDir(cCurrentPath, sizeof(cCurrentPath)))
	{
		std::cout << "Something went wrong when trying to obtain the current working directory" << std::endl;
		return;
	}

	cCurrentPath[sizeof(cCurrentPath) - 1] = '\0';

	printf ("The current working directory is %s \n", cCurrentPath);
	return;
}

