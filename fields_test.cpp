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

int main() 
{
	bool a = grid_test(); 
	
	return 0; 
	
}
