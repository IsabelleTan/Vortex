//
// Created by Isabelle Tan on 04-08-16.
//

#include <iostream>
#include <cmath>
#include "ADIdiffusion_test.h"
#include "ADIdiffusion.h"
#include "datapoints.h"


using namespace std;
const value_type EPSILON = 0.000001;

// A test function for the ThomasAlg() function that solves a linear system with a tridiagonal matrix
// IMPORTANT: Assuming that the k*n'th elements of the sub- and superdiagonal are 0.
bool ThomasAlg_test(){
    // Initialize parameters
    bool result = true;
    int N = 4;
    value_type * x = new value_type[N];
    value_type * control = new value_type[N]{1, 2, 3, 4};
    value_type * r = new value_type[N]{0,3,2,5};
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

    delete[] x;
    delete[] control;
    delete[] r;

    return result;
}

// A function to test the changeOrdering() function
bool changeOrdering_test(){
    bool result = true;

    // Allocate the arrays
    int N = 9;
    value_type * q_in = new value_type[N]{1,2,3,1,2,3,1,2,3};
    value_type * q_out = new value_type[N];
    value_type * control = new value_type[N]{1,1,1,2,2,2,3,3,3};

    // Perform the reordering
    changeOrdering(N, q_in, q_out);

    // Test the result
    for (int j = 0; j < N; ++j) {
        if(::abs(q_out[j] - control[j]) < EPSILON){
            cout << "q_out[j] = " << q_out[j] << " == " << control[j] << " = control[j]" << endl;
        } else {
            cout << "q_out[j] = " << q_out[j] << " != " << control[j] << " = control[j]" << endl;
            result = false;
        }
    }

    // Print the result
    if(result){
        cout << "Test succeeded!" << endl;
    } else {
        cout << "Test failed" << endl;
    }

    // Free memory
    delete[]  q_in;
    delete[] q_out;
    delete[] control;

    return result;
}

// A test function for the diffusion that generates 100 output files which we can animate
void ADI_test_output(){
    // Set parameters
    const value_type t_0 = 0;
    const value_type t_end = 10;
    const value_type dt = 0.1;
    int iter = (int)(t_end - t_0)/dt;
    const value_type dx = 0.1;
    const value_type v = 1;
    const int N = 2500;
    const int n = sqrt(N);

    // Allocate arrays
    value_type * q_0 = new value_type[N];
    value_type * q_new = new value_type[N];

    // Fill initial array (loop row wise)
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if(::abs(i-n/2)<10 && ::abs(j-n/2) < 10){
                q_0[i*n+j] = 1;
            } else{
                q_0[i*n+j] = 0;
            }
        }
    }

    // Assign a string for filename and write to file
    string filename = to_string(0) + ".txt";
    write_to_file(filename.c_str(), N, q_0);

    // Time loop
    for (int i = 1; i < iter; ++i) {
        // Compute one diffusion step
        ADI(N, q_0, q_new, dt, dx, v);

        // Assign a string for filename and write to file
        string filename = to_string(i) + ".txt";
        write_to_file(filename.c_str(), N, q_new);

        // Swap pointers q_0 to q_new
        value_type * temp = q_0;
        q_0 = q_new;
        q_new = temp;
    }

    // Free memory
    delete[] q_0;
    delete[] q_new;

    return;
}
