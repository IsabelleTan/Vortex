#include <iostream> 
#include <iomanip>  	
#include <string>
#include <chrono>

#include "morton.h"
#include "datapoints.h"

const double samplesize = 100.;
const unsigned int numberdatapoints = 1000000;

// mean of 100 randomly generated samples of size 10^6
double time_extent(const int N, const double* const x, const double* const y, double& xmin, double& ymin, double& ext){
	
	auto start = std::chrono::high_resolution_clock::now();
	auto end   = std::chrono::high_resolution_clock::now();
	double ext_time(0);
	
	for(size_t i(0); i<samplesize; i++){
		start = std::chrono::high_resolution_clock::now();
		extent(N, x, y, xmin, ymin, ext); 
		end = std::chrono::high_resolution_clock::now();
		ext_time += static_cast<std::chrono::duration<double>>(end-start).count();
	}
	
	return ext_time/samplesize;
}

double time_morton(const int N, const double* const x, const double* const y, const double xmin, const double ymin, const double ext, int* index){

	auto start = std::chrono::high_resolution_clock::now();
	auto end   = std::chrono::high_resolution_clock::now();
	double morton_time(0);
	
	for(size_t i=0; i<samplesize; i++){
		start = std::chrono::high_resolution_clock::now();
		morton(N, x, y, xmin, ymin, ext, index); 
		end = std::chrono::high_resolution_clock::now();
		morton_time += static_cast<std::chrono::duration<double>>(end-start).count();
	}
	
	return morton_time/samplesize;
}

double time_sort(const int N, int* index, int* keys){

	auto start = std::chrono::high_resolution_clock::now();
	auto end   = std::chrono::high_resolution_clock::now();
	double sort_time(0);
	
	for(size_t i=0; i<samplesize; i++){
		start = std::chrono::high_resolution_clock::now();
		sort(N, index, keys); 
		end = std::chrono::high_resolution_clock::now();
		sort_time += static_cast<std::chrono::duration<double>>(end-start).count();	
	}
	
	return sort_time/samplesize;
}

double time_reorder(const int N, const int* const keys, const double* const x, const double* const y, const double* const q, double* xsorted, double* ysorted, double* qsorted){
	
	auto start = std::chrono::high_resolution_clock::now();
	auto end   = std::chrono::high_resolution_clock::now();
	double reorder_time(0);
	
	for(size_t i=0; i<100; i++){
		start = std::chrono::high_resolution_clock::now();
		reorder(N, keys, x, y, q, xsorted, ysorted, qsorted); 
		end = std::chrono::high_resolution_clock::now();
		reorder_time += static_cast<std::chrono::duration<double>>(end-start).count();
	}
	
	return reorder_time/samplesize;
}

void time_all(){														// time all functions and display times 
	// declare the variables and vectors needed
	double *x, *y, *q, *xsorted, *ysorted, *qsorted;
	int *index, *keys; 
	double xmin, ymin, ext; 
	std::chrono::time_point< std::chrono::high_resolution_clock > start , end;
		
	// fill x and y up and memory-align them as well 
	load_data_random(numberdatapoints, &x, &y);
	load_data_random(numberdatapoints, &q);
	
	// memory-align the other vectors
	posix_memalign( (void**)&xsorted, 32, sizeof(double)*numberdatapoints);
	posix_memalign( (void**)&ysorted, 32, sizeof(double)*numberdatapoints);
	posix_memalign( (void**)&qsorted, 32, sizeof(double)*numberdatapoints);
	posix_memalign( (void**)&index  , 32, sizeof(int)  *numberdatapoints);
	posix_memalign( (void**)&keys   , 32, sizeof(int)  *numberdatapoints);
	
	// fill up the "keys"-vector
	for(unsigned int i(0); i<numberdatapoints; i++)
		keys[i] = i; 
		
	// execute functions and perform time measurements
	const double ext_time(time_extent(numberdatapoints, x, y, xmin, ymin, ext));	
	const double morton_time(time_morton(numberdatapoints, x, y, xmin, ymin, ext, index));	
	const double sort_time(time_sort(numberdatapoints, index, keys));
	const double reorder_time(time_reorder(numberdatapoints, keys, x, y, q, xsorted, ysorted, qsorted));
	const double total_time = ext_time + morton_time + sort_time + reorder_time ;
	
	// display time measurements and percentage of total time
	std::cout << "For " << numberdatapoints << " data points averaged over " << samplesize << " samples:" << std::endl
			  << std::setw(25) << "function extent (4t) : "  << ext_time*1000     << " ms     = " 
														<< (ext_time/total_time)*100. << "% of total time" 
														<< std::endl
			  << std::setw(25) << "function morton : "  << morton_time*1000  << " ms     = " 
														<< (morton_time/total_time)*100. << "% of total time" 
														<< std::endl
			  << std::setw(25) << "function sort : "    << sort_time*1000    << " ms     = " 
														<< (sort_time/total_time)*100. << "% of total time" 
														<< std::endl
			  << std::setw(25) << "function reorder : " << reorder_time*1000 << " ms     = " 
														<< (reorder_time/total_time)*100. << "% of total time" 
														<< std::endl << std::endl; 
	
	// liberate the allocated memory  	
	free(x); 
	free(y);
	free(q);
	free(xsorted); 
	free(ysorted); 
	free(index); 
	free(keys);	
}

int main(){
	
	time_all();
	
	return 0;
	
}
