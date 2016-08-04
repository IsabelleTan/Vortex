//
// Created by Isabelle Tan on 04-08-16.
//

#include <iostream>
#include "ADIdiffusion_test.h"
#include "ADIdiffusion.h"


using namespace std;
const value_type EPSILON = 0.000001;

// A test function for the ThomasAlg() function that solves a linear system with a tridiagonal matrix
bool ThomasAlg_test(){
    // Initialize parameters
    bool result = true;
    int N = 4;
    value_type * x = new value_type[N];
    value_type * control = new value_type[N]{1, 2, 3, 4};
    value_type * r = new value_type[N]{0,0,0,5};
    const value_type a = -1;
    const value_type b = 2;
    const value_type c = -1;

    // Perform the computation
    ThomasAlg(N, a, b, c, x, r);

    // Test the result
    for (int j = 0; j < N; ++j) {
        if(::abs(x[j] - control[j]) < EPSILON){
            cout << "x[j] = " << x[j] << " == " << control[j] << " = control[j]" << endl;
        } else {
            cout << "x[j] = " << x[j] << " != " << control[j] << " = control[j]" << endl;
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
