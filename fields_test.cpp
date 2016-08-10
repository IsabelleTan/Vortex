//
// Created by Isabelle Tan on 28-07-16.
//

#include "fields.h"
#include <iostream>

// A test function for the grid() function
bool grid_test(){
    bool result = true;

    // Set parameters
    const int N = 9;
    const value_type xRange = 8;
    const value_type yRange = 8;

    // Create the test arrays
    const value_type control_x[N] = {-4, 0, 4, -4, 0, 4, -4, 0, 4};
    const value_type control_y[N] = {4, 4, 4, 0, 0, 0, -4, -4, -4};
    value_type x[N];
    value_type y[N];

    // Compute
    grid(N, x, y, xRange, yRange);

    // Check
    for (int i = 0; i < N; ++i) {
        if (x[i]!=control_x[i] || y[i]!=control_y[i]){
            result = false;
            std::cout << " for i = " << i << " x[i] = " << x[i] << " != " << control_x[i] 
					  << " = control_x[i] \n or \n y[i] = " << y[i] << " != " << control_y[i] 
					  << " = control_y[i]" << std::endl;
        }
    }


    // Print the result
    if(result){
        std::cout << "Test succeeded!" << std::endl;
    } else {
        std::cout << "Test failed" << std::endl;
    }


    return result;
}

bool arrayToMatrix_test(){
    bool result = true;
    int N = 100;
    int n = sqrt(N);

    // Prepare arguments
    MatrixXd matrix(n,n);
    value_type * array = new value_type[100];

    // Fill array with values
    for (int i = 0; i < 100; ++i) {
        array[i] = i;
    }

    // Convert to matrix
    arrayToMatrix(array, matrix, false);

    // Print content of matrix
    for (int j = 0; j < sqrt(N); ++j) {
        for (int i = 0; i < sqrt(N); ++i) {
            std::cout << "Array[" << j << "," << i << "] = " << array[j*(int)sqrt(N) + i] << std::endl;
            std::cout << "Matrix[" << j << "," << i << "] = " << matrix(j, i) << std::endl;
        }
    }

    delete[] array;
    return result;
}


int main() 
{
	arrayToMatrix_test();
	
	return 0; 
	
}
