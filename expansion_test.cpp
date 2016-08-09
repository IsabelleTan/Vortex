
#include "morton.h"
#include "expansion.h"
#include <iostream>

value_type epsilon = 0.000001;

/*
 * A function to test and time the p2e kernel
 */

value_type timep2e(int N, int order){

    // Create arrays with particles
    value_type *x = new value_type[N];
    value_type *y = new value_type[N];
    value_type *mass = new value_type[N];

    value_type xCom = 0.5;
    value_type yCom = 0.5;

    initialize(N, x, y, mass);

    value_type *expansion = new value_type[2*order]{0};

	std::chrono::time_point< std::chrono::high_resolution_clock > start , stop;
	start = std::chrono::high_resolution_clock::now();					// start time 
    
    p2e(x, y, mass, N, order, xCom, yCom, expansion, expansion + order);
    
    stop = std::chrono::high_resolution_clock::now();        			// record end time
	const double time = static_cast<std::chrono::duration<double>>(stop-start).count();

    delete[] x;
    delete[] y;
    delete[] mass;
    delete[] expansion;

    return time;
}

/*
 * A function to test the p2e kernel
 */
bool testp2e(){
    bool result = true;
    int order = 2;
    int N = 4;

    value_type x[4] = {1, 2, 3, 4};
    value_type y[4] = {4, 3, 2, 1};
    value_type mass[4] = {5, 6, 7, 8};

    value_type xCom= 2.5;
    value_type yCom = 2.5;

    value_type expansion[4] = {0};

    p2e(x, y, mass, N, order, xCom, yCom, expansion, expansion + order);

    value_type control[4] = {-70, -25, -60, -130};

    for (int i = 0; i < 4; ++i) {
        if (abs(expansion[i] - control[i]) > epsilon){
            std::cout << "expansion[" << i << "] = " << expansion[i] << " != control[" << i << "] = " << control[i] << std::endl;
            result = false;
        } else {

            std::cout << "expansion[" << i << "] = " << expansion[i] << " = control[" << i << "] = " << control[i] << std::endl;
        }
    }

    if (result){
        std::cout << "Test succeeded" << std::endl;
    }

    return result;
}

/*
 * A function to test and time the e2p kernel
 */
value_type timee2p(int N, int order){

    // Create arrays with particles
    value_type *x = new value_type[N];
    value_type *y = new value_type[N];
    value_type *mass = new value_type[N];

    value_type xCom = 0.5;
    value_type yCom = 0.5;

    initialize(N, x, y, mass);

    value_type *expansion = new value_type[2*order]{0};

    p2e(x, y, mass, N, order, xCom, yCom, expansion, expansion + order);

    //Prepare arguments
    value_type xtarget= 3;
    value_type ytarget = 1;
    value_type q = 0;

    // Compute g (sum of masses)
    for (int i = 0; i < N; ++i) {
        q+= mass[i];
    }

    std::chrono::time_point< std::chrono::high_resolution_clock > start, stop;
	start = std::chrono::high_resolution_clock::now();		// start time 
    
    value_type streamfunction = e2p(xtarget, ytarget, q, order, expansion, expansion + order);
    
    stop = std::chrono::high_resolution_clock::now();        	// record end time
	const double time = static_cast<std::chrono::duration<double>>(stop-start).count();

    delete[] x;
    delete[] y;
    delete[] mass;
    delete[] expansion;

    return time;
}

/*
 * A function to test the e2p kernel
 */
bool teste2p(){
    bool result = true;
    int order = 2;
    int N = 4;

    value_type x[4] = {1, 2, 3, 4};
    value_type y[4] = {4, 3, 2, 1};
    value_type mass[4] = {5, 6, 7, 8};

    value_type xCom= 2.5;
    value_type yCom = 2.5;

    value_type x_target = 3;
    value_type y_target = 1;

    value_type expansion[2*order] = {0};

    p2e(x, y, mass, N, order, xCom, yCom, expansion, expansion + order);

    value_type streamfunction = e2p(x_target, y_target, 26, order, expansion, expansion + order);

    value_type control = -6.8663938;

    if (abs(streamfunction - control) > epsilon){
        std::cout << "streamfunction = " << streamfunction << " != control = " << control << std::endl;
        result = false;
    } else {

        std::cout << "streamfunction = " << streamfunction << " = control = " << control << std::endl;
    }


    if (result){
        std::cout << "Test succeeded" << std::endl;
    }

    return result;
}

/*
 * A function for the convergence analysis of the p2e+e2p kernels
 */
void convergenceAnalysis(int p_max){
    int N = 100000;

    // Create arrays with particles
    value_type *x = new value_type[N];
    value_type *y = new value_type[N];
    value_type *mass = new value_type[N];

    initialize(N, x, y, mass);

    //Prepare arguments
    value_type xtarget= 3;
    value_type ytarget = 1;
    value_type q = 0;
    value_type xCom = 0;
    value_type yCom = 0;

    // Compute q (sum of masses) and xCom and yCom (center of mass)
    for (int i = 0; i < N; ++i) {
        q+= mass[i];
        xCom+= mass[i]*x[i];
        yCom+= mass[i]*y[i];
    }
    xCom/=q;
    yCom/=q;


    // Compute the streamfunction for different orders and print the resulting value
    for (int p = 3; p <= p_max; ++p) {
        value_type *expansion = new value_type[2*p]{0};

        p2e(x, y, mass, N, p, xCom, yCom, expansion, expansion + p);

        std::cout << e2p(xtarget, ytarget, q, p, expansion, expansion + p) << std::endl;

        delete[] expansion;
    }


    delete[] x;
    delete[] y;
    delete[] mass;
}

/*
 * A function to time the p2p kernel
 */
value_type timep2p(int N){

    // Create arrays with particles
    value_type *x = new value_type[N];
    value_type *y = new value_type[N];
    value_type *mass = new value_type[N];

    value_type xtarget = 3;
    value_type ytarget = 1;

    initialize(N, x, y, mass);

	std::chrono::time_point< std::chrono::high_resolution_clock > start , stop;
	start = std::chrono::high_resolution_clock::now();		// start time 

    value_type streamfunction = p2p(x, y, mass, N, xtarget, ytarget);

	stop = std::chrono::high_resolution_clock::now();        	// record end time
	const double time = static_cast<std::chrono::duration<double>>(stop-start).count();

    delete[] x;
    delete[] y;
    delete[] mass;

    return time;

}

/*
 * A function to test the p2p kernel
 */
bool testp2p(){
    bool result = true;
    int order = 2;
    int N = 4;

    value_type x[4] = {1, 2, 3, 4};
    value_type y[4] = {4, 3, 2, 1};
    value_type mass[4] = {5, 6, 7, 8};

    value_type xCom= 2.5;
    value_type yCom = 2.5;

    value_type x_target = 3;
    value_type y_target = 1;


    value_type streamfunction = p2p(x, y, mass, N, x_target, y_target);

    value_type control = 11.240687;

    if (abs(streamfunction - control) > epsilon){
        std::cout << "streamfunction = " << streamfunction << " != control = " << control << std::endl;
        result = false;
    } else {

        std::cout << "streamfunction = " << streamfunction << " = control = " << control << std::endl;
    }


    if (result){
        std::cout << "Test succeeded" << std::endl;
    }

    return result;
}

/*
 * A function to time something N times and compute the average and variance
 */
void time(int N){
    value_type *times = new value_type[N];

    // Track the time N times
    for (int i = 0; i < N; ++i) {
        times[i] = timee2p(100000,20);
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

int main() {
	
	value_type a = timep2e(1000000, exp_order);
	bool b = testp2e();
	a = timee2p(1000000, exp_order);
	b = teste2p();
	convergenceAnalysis(exp_order);
	a = timep2p(1000000);
	b = testp2p();
	
	return 0;
	
}


