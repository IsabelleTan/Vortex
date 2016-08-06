//
// Created by Isabelle Tan on 27-07-16.
//

#include "solver.h"
#include "solver_test.h"
#include <iostream>

using namespace std;

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
        control_u[i] = -4; //TODO this still needs to be computed to be tested
        control_v[i] = -1;
    }

    // Perform computation
    velocity(N, 1, u, v, phi);

    // Test the result
    for (int j = 0; j < N; ++j) {
        if(u[j] == control_u[j] && v[j] == control_v[j]){
            cout << "u[j] = " << u[j] << " == " << control_u[j] << " = control_u[j]" << " and" << endl;
            cout << "v[j] = " << v[j] << " == " << control_v[j] << " = control_v[j]\n" << endl;
        } else {
            cout << "u[j] = " << u[j] << " != " << control_u[j] << " = control_u[j]" << " or" << endl;
            cout << "v[j] = " << v[j] << " != " << control_v[j] << " = control_v[j]\n" << endl;
            result = false;
        }
    }

    // Print the result
    if(result){
        cout << "Test succeeded!" << endl;
    } else {
        cout << "Test failed" << endl;
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
            cout << "q[j] = " << q[j] << " == " << control_q[j] << " = control_q[j]" << endl;
        } else {
            cout << "q[j] = " << q[j] << " != " << control_q[j] << " = control_q[j]" << endl;
            result = false;
        }
    }

    // Print the result
    if(result){
        cout << "Test succeeded!" << endl;
    } else {
        cout << "Test failed" << endl;
    }

    return result;
}