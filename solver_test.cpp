//
// Created by Isabelle Tan on 27-07-16.
//

#include "solver.h"
#include "datapoints.h"
#include <iostream>
#include <chrono>
#include <cmath>

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

    delete[] control_u;
    delete[] control_v;
    return result;
}

/*
 * A function to time the velocity.
 */
value_type velocity_time(int N){
    value_type time;
    value_type * phi;
    value_type * u;
    value_type * v;

    posix_memalign( (void**) &phi, 32, sizeof(value_type)*N );
    posix_memalign( (void**) &u, 32, sizeof(value_type)*N );
    posix_memalign( (void**) &v, 32, sizeof(value_type)*N );

    // Prepare time variables
    auto start = std::chrono::high_resolution_clock::now();
    auto end   = std::chrono::high_resolution_clock::now();

    // Initialize phi with random values
    load_data_random(N, &phi);

    start = std::chrono::high_resolution_clock::now();
    velocity(N, 1, u, v, phi);
    end = std::chrono::high_resolution_clock::now();
    time = static_cast<std::chrono::duration<double>>(end-start).count();


    delete[] phi;
    delete[] u;
    delete[] v;

    return time;
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

    delete[] u;
    delete[] v;
    delete[] control_q;
    return result;
}

/*
 * A function to time the vorticity.
 */
value_type vorticity_time(int N){
    value_type time;
    value_type * q = new value_type[N];
    value_type * u = new value_type[N];
    value_type * v = new value_type[N];

    posix_memalign( (void**) &q, 32, sizeof(value_type)*N );
    posix_memalign( (void**) &u, 32, sizeof(value_type)*N );
    posix_memalign( (void**) &v, 32, sizeof(value_type)*N );

    // Prepare time variables
    auto start = std::chrono::high_resolution_clock::now();
    auto end   = std::chrono::high_resolution_clock::now();

    // Initialize phi with random values
    load_data_random(N, &u, &v);

    start = std::chrono::high_resolution_clock::now();
    vorticity(N, 1, u, v, q);
    end = std::chrono::high_resolution_clock::now();
    time = static_cast<std::chrono::duration<double>>(end-start).count();


    delete[] q;
    delete[] u;
    delete[] v;

    return time;
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

    delete[] u;
    delete[] v;
    delete[] x;
    delete[] y;
    delete[] control_x;
    delete[] control_y;
    return result;
}

/*
 * A function to time the advection.
 */
value_type advection_time(int N){
    value_type time;
    value_type * x = new value_type[N];
    value_type * y = new value_type[N];
    value_type * u = new value_type[N];
    value_type * v = new value_type[N];


    posix_memalign( (void**) &x, 32, sizeof(value_type)*N );
    posix_memalign( (void**) &y, 32, sizeof(value_type)*N );
    posix_memalign( (void**) &u, 32, sizeof(value_type)*N );
    posix_memalign( (void**) &v, 32, sizeof(value_type)*N );

    // Prepare time variables
    auto start = std::chrono::high_resolution_clock::now();
    auto end   = std::chrono::high_resolution_clock::now();

    // Initialize phi with random values
    load_data_random(N, &x, &y);
    load_data_random(N, &u, &v);

    start = std::chrono::high_resolution_clock::now();
    advection(N, 1, u, v, x, y);
    end = std::chrono::high_resolution_clock::now();
    time = static_cast<std::chrono::duration<double>>(end-start).count();


    delete[] x;
    delete[] y;
    delete[] u;
    delete[] v;

    return time;
}


/*
 * A function to time something N times and compute the average and variance
 */
void time(int N){
    value_type *times = new value_type[N];

    // Track the time N times
    for (int i = 0; i < N; ++i) {
        times[i] = velocity_time(1000000);
    }

    // Compute the average
    value_type mean = 0;
    for (int j = 0; j < N; ++j) {
        mean += times[j];
    }
    mean/=N;

    // Compute the variance
    value_type var = 0;
    for (int k = 0; k < N; ++k) {
        var += pow(times[k] - mean, 2);
    }
    var/=N;

    // Print the results
    std::cout << "Timed " << N << " times." << std::endl;
    std::cout << "Mean = " << mean << std::endl;
    std::cout << "Variance = " << var << " is " << var/mean * 100 << "% of the mean" << std::endl;

    delete[] times;
}

int main()
{
	
	//bool a = velocity_test();
	
	//bool b = vorticity_test();

    //bool c = advection_test();

    time(100);
	
	return 0; 
	
} 
