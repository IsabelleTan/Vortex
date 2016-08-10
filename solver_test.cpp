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
    value_type * control_u = new value_type[N]{0, 0, 0, 0, 0, -4, -4, 0, 0, -4, -4, 0, 0, 0, 0, 0};
    value_type * control_v = new value_type[N]{0, 0, 0, 0, 0, -1, -1, 0, 0, -1, -1, 0, 0, 0, 0, 0};

    // Initialize
    for (int i = 0; i < N; ++i) {
        phi[i] = i;
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
    int N = 9;

    // Allocate and initialize arrays
    value_type q[N];
    value_type * u = new value_type[N]{2,3,2,3,2,3,2,3,2};
    value_type * v = new value_type[N]{1,2,3,4,5,6,7,8,9};
    value_type * control_q = new value_type[N]{-0.5, 0, -2.5, 2.5, 1, -2.5, 5.5, 2, -2.5};

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

/*
 * A function to test the advection() function.
 */
bool advection_test(){
    // Set parameters
    bool result = true;
    int N = 9;

    // Allocate and initialize arrays
    value_type * u = new value_type[N]{2,4,2,4,2,4,2,4,2};
    value_type * v = new value_type[N]{8,6,8,6,8,6,8,6,8};
    value_type * x = new value_type[N]{0,0,0,0,0,0,0,0,0};
    value_type * y = new value_type[N]{0,0,0,0,0,0,0,0,0};
    value_type * control_x = new value_type[N]{1,2,1,2,1,2,1,2,1};
    value_type * control_y = new value_type[N]{4,3,4,3,4,3,4,3,4};

    // Perform computation
    advection(N, 0.5, u, v, x, y);

    // Test the result
    for (int j = 0; j < N; ++j) {
        if(x[j] == control_x[j] && y[j] == control_y[j]){
            std::cout << "x[j] = " << x[j] << " == " << control_x[j] << " = control_x[j]" << " and" << std::endl;
            std::cout << "y[j] = " << y[j] << " == " << control_y[j] << " = control_y[j]\n" << std::endl;
        } else {
            std::cout << "x[j] = " << x[j] << " != " << control_x[j] << " = control_x[j]" << " or" << std::endl;
            std::cout << "y[j] = " << y[j] << " != " << control_y[j] << " = control_y[j]\n" << std::endl;
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
	
	//bool a = velocity_test();
	
	//bool b = vorticity_test();

    //bool c = advection_test();
	
	return 0; 
	
} 
