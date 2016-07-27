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
        control_u[i] = -4;
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
            cout << "or" << endl;
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