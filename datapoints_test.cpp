#include "datapoints.h"
#include <string>
#include <iostream>
#include <cmath>

void testLoadingAndWriting(){
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
}


int main(){
	// Create coordinate vectors
	int N = 10000;
	value_type * x;
	value_type * y;
	value_type * q = new value_type[N];

	// Fill with random data
	load_data_random(N, &x, &y);

	// Define filenames
	std::string filenameX = "0_X.txt";
	std::string filenameY = "0_Y.txt";
	std::string filenameQ = "0_Q.txt";

	// Compute some multivariate function
	for (int t = 0; t < 10; ++t) {
		for (int i = 0; i < N; ++i) {
			q[i] = sin(2 * x[i] + t);
			std::cout << "q[i] = " << q[i] << std::endl;

		}
		filenameQ = std::to_string(t) + "_Q.txt";
		filenameX = std::to_string(t) + "_X.txt";
		filenameY = std::to_string(t) + "_Y.txt";

		// Write to files
		write_to_file(filenameQ.c_str(), N, q);
		write_to_file(filenameX.c_str(), N, x);
		write_to_file(filenameY.c_str(), N, y);
	}







	delete[] q;
	return 0;
}
