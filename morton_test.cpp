#include <iostream> 
#include <iomanip>  	
#include <bitset>
#include <chrono>

#include "morton.h"
#include "datapoints.h"

void test_extent()
{
	const int N = 6;
	double xx[N] = {89, 83, 13, 0.23, 84, 5.66} ;
	double yy[N] = {45.4, 5.4, 67.9, 4.66, 43.1, 23.4};
	double x_min, y_min, ext;
	
	extent(N, xx, yy, x_min, y_min, ext);
	
	std::cout << "Extent:" << std::endl
			  << "xmin = " << x_min << std::endl
			  << "ymin = " << y_min << std::endl
			  << "ext = " << ext << std::endl;
}

void test_morton()
{
	
	const unsigned int depth(15);	
	const int N = 6;
	double xx[N];
	double yy[N];
	unsigned int ind[N];
	
	xx[0] = 0.; 	yy[0] = 0.;
	xx[1] = 100.; 	yy[1] = 100.;
	xx[2] = 0.9;	yy[2] = 0.3;
	xx[3] = 98.2; 	yy[3] = -99.1;
	xx[4] = 1.; 	yy[4] = 99.1;
	xx[5] = 96.7; 	yy[5] = 0.3;	
	
		
	double x_min, y_min, ext;	
	extent(N, xx, yy, x_min, y_min, ext);
	morton(N, xx, yy, x_min, y_min, ext, ind, depth);

	std::cout << " extent = " << ext << std::endl;
	std::cout << " Morton indices: " << std::endl;
	for(size_t i(0); i<N; ++i) {
		std::cout << "[" << i << "] (" << std::setw(4) << xx[i]
		<< ", " << std::setw(4) << yy[i]
		<< std::setw(12) << ") <-----> " << (std::bitset<30>(ind[i])).to_string()
		<< std::endl;
	}
	for (int j = 0; j < N; ++j) {
		std::cout << "ind[" << j << "] = " << ind[j] << std::endl;
	}
} 

void test_sort()
{
	const unsigned int depth(15);
	const int N = 6;
	double xx[N];
	double yy[N];
	unsigned int ind[N]; // TODO changed this from int to unsigned int
	unsigned int key[N];
	
	xx[0] = 0.; 	yy[0] = 0.;
	xx[1] = 100.; 	yy[1] = 100.;
	xx[2] = 30.4;	yy[2] = 48;
	xx[3] = 4.2; 	yy[3] = 67.1;
	xx[4] = 85.45; 	yy[4] = 32.1;
	xx[5] = 33.7; 	yy[5] = 16.3;	
	
	double x_min, y_min, ext;	
	extent(N, xx, yy, x_min, y_min, ext);
	morton(N, xx, yy, x_min, y_min, ext, ind, depth);


	std::cout << " extent = " << ext << std::endl;
	std::cout << " Morton indices: " << std::endl;
	for(size_t i(0); i<N; ++i){
		std::cout << "[" << i << "] (" << std::setw(4) << xx[i] 
				  << ", " << std::setw(4) << yy[i] 
		          << std::setw(12) << ") <-----> " << (std::bitset<30>(ind[i])).to_string() 
		          << std::endl;
	}	
	
	for(size_t i(0); i<N; i++){
		key[i] = i;
	}
	
	sort(N, ind, key);
	std::cout << "now let's sort !" << std::endl; 
	for(size_t i(0); i<N; ++i){
		std::cout << (std::bitset<30>(ind[i])).to_string() 
				  << ",    " << ind[i] 
				  << ",    " << key[i] 
		          << std::endl;
	}	
}

void test_all(int n)
{
	// declare the variables and vectors needed
	const int depth(15);
	double *x, *y, *q, *xsorted, *ysorted, *qsorted;
	unsigned int *index;
	unsigned int *keys;
	double xmin, ymin, ext; 
	std::chrono::time_point< std::chrono::high_resolution_clock > start , end;
		
	// fill x and y up and memory-align them as well 
	if(n<0){													// if n<0, then: read data from a file 
		load_data_from_file("test_datasets/test-1.blob", n, &x, &y);
		load_data_random(n, &q);
	} else {													// if n>=0, then: set random values
		load_data_random(n, &x, &y);
		load_data_random(n, &q);
	}
		
	// memory-align the other vectors
	posix_memalign( (void**)&xsorted, 32, sizeof(double)*n);
	posix_memalign( (void**)&ysorted, 32, sizeof(double)*n);
	posix_memalign( (void**)&qsorted, 32, sizeof(double)*n);
	posix_memalign( (void**)&index  , 32, sizeof(int)  *n);
	posix_memalign( (void**)&keys   , 32, sizeof(int)  *n);
	
	// fill up the "keys"-vector
	for(unsigned int i(0); i<n; i++)
		keys[i] = i; 
	
	// execute functions and perform time measurements
	start = std::chrono::high_resolution_clock::now();
	extent(n, x, y, xmin, ymin, ext); 
	end = std::chrono::high_resolution_clock::now();
	const double ext_time = static_cast<std::chrono::duration<double>>(end-start).count();

	start = std::chrono::high_resolution_clock::now();
	morton(n, x, y, xmin, ymin, ext, index, depth); 
	end = std::chrono::high_resolution_clock::now();
	const double morton_time = static_cast<std::chrono::duration<double>>(end-start).count();
	
	start = std::chrono::high_resolution_clock::now();
	sort(n, index, keys); 
	end = std::chrono::high_resolution_clock::now();
	const double sort_time = static_cast<std::chrono::duration<double>>(end-start).count();
	
	start = std::chrono::high_resolution_clock::now();
	reorder(n, keys, x, y, q, xsorted, ysorted, qsorted); 
	end = std::chrono::high_resolution_clock::now();
	const double reorder_time = static_cast<std::chrono::duration<double>>(end-start).count();
	
	const double total_time = ext_time + morton_time + sort_time + reorder_time ;
	
	// display time measurements and percentage of total time
	std::cout << "For " << n << " data points:" << std::endl
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

int main() 
{
	// --> testing individual functions
	//test_morton();
	//test_extent();
	//test_sort();
	
	// --> perform check of all functions
	//test_all(10000000);
	//test_all(-1);
	
	//! results (serial)
	/*
	For 10000000 data points (random):
       function extent : 182.953 ms     = 3.37604% of total time
       function morton : 223.03 ms     = 4.11559% of total time
         function sort : 4707.28 ms     = 86.8639% of total time
      function reorder : 305.88 ms     = 5.64442% of total time

	For 9874432 data points (file "test-2.blob"):
       function extent : 106.606 ms     = 2.15128% of total time
       function morton : 219.766 ms     = 4.43482% of total time
         function sort : 4581.12 ms     = 92.4459% of total time
      function reorder : 47.9701 ms     = 0.968025% of total time
     */
     
     //! results (parallel)
     /*
	 For 10000000 data points (random):
   function extent (p) : 93.8466 ms     = 3.93975% of total time
       function morton : 102.352 ms     = 4.29679% of total time
         function sort : 2016.68 ms     = 84.6616% of total time
      function reorder : 169.169 ms     = 7.10185% of total time

	 For 9874432 data points (file "test-2.blob"):
   function extent (p) : 56.6856 ms     = 2.28677% of total time
       function morton : 78.9418 ms     = 3.18462% of total time
         function sort : 2318.47 ms     = 93.5302% of total time
      function reorder : 24.7488 ms     = 0.9984% of total time
     */
     
     //! Analysis
     //! function "sort" takes by far the most time, that's the one we should work on

		
	return 0; 
}


