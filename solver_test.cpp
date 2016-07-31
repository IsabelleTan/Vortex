//
// Created by Isabelle Tan on 27-07-16.
//

#include "solver.h"
#include <iostream>

/*
 * A function to test the velocity() function
 */
bool velocity_test(){
    // Set parameters
    bool result = true;
    int N = 16;

    // Allocate arrays
    value_type phi[N];
    value_type u[N];
    value_type v[N];
    value_type control_u[N];
    value_type control_v[N];

    // Initialize
    for (int i = 0; i < N; ++i) {
        phi[i] = i;
        control_u[i] = -4;
        control_v[i] = -1;
    }

    // Perform computation
    velocity(N, 1, u, v, phi);

    // Test the result
    for (int j = 0; j < N; ++j) {
        if(u[j] == control_u[j] && v[j] == control_v[j]){
            std::cout << "u[j] = " << u[j] << " == " << control_u[j] << " = control_u[j]" << " and" << std::endl;
            std::cout << "v[j] = " << v[j] << " == " << control_v[j] << " = control_v[j]\n" << std::endl;
        } else {
            std::cout << "u[j] = " << u[j] << " != " << control_u[j] << " = control_u[j]" << " or" << std::endl;
            std::cout << "v[j] = " << v[j] << " != " << control_v[j] << " = control_v[j]\n" << std::endl;
            result = false;
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

/*
 * A function to test the vorticity() function
 */
bool vorticity_test(){
    // Set parameters
    bool result = true;
    int N = 16;

    // Allocate arrays
    value_type q[N];
    value_type u[N];
    value_type v[N];
    value_type control_q[N];

    // Initialize
    for (int i = 0; i < N; ++i) {
        u[i] = i;
        v[i] = i + 1;
        control_q[i] = 5;
    }

    // Perform computation
    vorticity(N, 1, u, v, q);

    // Test the result
    for (int j = 0; j < N; ++j) {
        if(q[j] == control_q[j]){
            std::cout << "q[j] = " << q[j] << " == " << control_q[j] << " = control_q[j]" << std::endl;
        } else {
            std::cout << "q[j] = " << q[j] << " != " << control_q[j] << " = control_q[j]" << std::endl;
            result = false;
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
	
	bool a = velocity_test(); 
	
	bool b = vorticity_test(); 
	
	return 0; 
	
} 
