//
// Created by Isabelle Tan on 04-08-16.
//

#include <iostream>
#include <cmath>
#include "ADIdiffusion_test.h"
#include "ADIdiffusion.h"
#include "datapoints.h"
#include "fields.h"

using namespace Eigen;
using namespace std;
const value_type EPSILON = 0.000001;

// A test function for the ThomasAlg() function that solves a linear system with a tridiagonal matrix
bool ThomasAlg_test(){
    // Initialize parameters
    bool result = true;
    int N = 4;
    MatrixXd x(4,1);
    MatrixXd control(4,1);
    control << 1, 2, 3, 4;
    MatrixXd r(4,1);
    r << 0,0,0,5;
    const value_type a = -1;
    const value_type b = 2;
    const value_type c = -1;

    // Perform the computation
    ThomasAlg(N, a, b, c, x, r);

    // Test the result
    for (int j = 0; j < N; ++j) {
        if(::abs(x(j,0) - control(j,0)) < EPSILON){
            cout << "x(j,0) = " << x(j,0) << " == " << control(j,0) << " = control(j,0)" << endl;
        } else {
            cout << "x(j,0) = " << x(j,0) << " != " << control(j,0) << " = control(j,0)" << endl;
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

// A test function for the diffusion that generates 100 output files which we can animate
void ADI_test_output(){
    // Set parameters
    const value_type t_0 = 0;
    const value_type t_end = 10;
    const value_type dt = 0.1;
    int iter = (int)(t_end - t_0)/dt;
    const value_type dx = 1;
    const value_type v = 1;
    const int N = 2500;
    const int n = sqrt(N);

    // Allocate array
    value_type * q_0 = new value_type[N];

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

    // Put the array in a matrix
    MatrixXd q_0_mat(n,n);
    arrayToMatrix(q_0, q_0_mat);
    MatrixXd q_new_mat(n,n);
    value_type * q_new = new value_type[N];

    // Assign a string for filename and write to file
    string filename = to_string(0) + ".txt";
    write_to_file(filename.c_str(), N, q_0);

    // Time loop
    for (int i = 1; i < iter; ++i) {
        // Compute one diffusion step
        ADI(q_0_mat, q_new_mat, dt, dx, v);
        matrixToArray(q_new, q_new_mat);

        // Assign a string for filename and write to file
        string filename = to_string(i) + ".txt";
        write_to_file(filename.c_str(), N, q_new);

        // Set q_0 to be q_new_mat for the next timestep
        q_0_mat = q_new_mat;
    }

    // Free memory
    delete[] q_0;
    delete[] q_new;

    return;
}
