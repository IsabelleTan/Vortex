#include <iostream>
#include <chrono>
#include "morton.h"
#include "quadtree.h"

// A function that times the build function and writes the output to a file named "times.txt".
value_type timeBuild(int N){
    // Prepare parameters
    int depth = 16;
    int k = 8;

    value_type *x = new value_type[N];
    value_type *y = new value_type[N];
    value_type *mass = new value_type[N];

    initialize(N, x, y, mass);

    value_type *xsorted = new value_type[N];
    value_type *ysorted = new value_type[N];
    value_type *mass_sorted = new value_type[N];

    // Allocate the tree array containing all the nodes
    int maxNodes = (int) std::min((float)8 * N / k, (float)pow(4, depth));
    Node* tree = new Node[maxNodes];

    std::cout << "Building the tree for N = " << N << std::endl;

	std::chrono::time_point< std::chrono::high_resolution_clock > start , stop; // declare variables required for timing
	start = std::chrono::high_resolution_clock::now();							// start time

    build(x, y, mass, N, k, xsorted, ysorted, mass_sorted, tree, depth);
    stop = std::chrono::high_resolution_clock::now();        					// record end time
	const double time = static_cast<std::chrono::duration<double>>(stop-start).count();
    std::cout << "N = " << N << " \tElapsed time: " << time << std::endl;

    delete[] tree;
    delete[] x;
    delete[] y;
    delete[] mass;

    delete[] xsorted;
    delete[] ysorted;
    delete[] mass_sorted;


    return time;
}

int main(){
	
	value_type a = timeBuild(1000000);
	
	return 0;
	
}
