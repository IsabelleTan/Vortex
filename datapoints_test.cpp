#include "datapoints.h"
#include <string>
#include <iostream>

int main(){
	
	int N;
	value_type *x, *y, *z;

	load_data_from_file("test_datasets/test-0.blob", N, &x, &y);
	
	std::cout << "Results: " << std::endl
			  << "N = " << N << std::endl
			  << "Values: " << N << std::endl;
	for(size_t i=0; i<15; i++){
		std::cout << x[i] << std::endl;
	}
	
	write_to_file("data-42.dat", N/2, x);
	
	load_data_from_file("test_datasets/test-0.blob", N, &z, &y);
	
	std::cout << "Results: " << std::endl
			  << "N = " << N << std::endl
			  << "Values: " << N << std::endl;
	for(size_t i=0; i<15; i++){
		std::cout << z[i] << std::endl;
	}

	return 0;
}
