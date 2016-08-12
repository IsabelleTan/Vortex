//
// Created by Isabelle Tan on 27-07-16.
//

#include <iostream>
#include "fields.h"
#include "datapoints.h"
#include "simulation.h"

/*
 * This function creates a regular 2D grid of N (x,y) coordinates with domain (xRange,yRange) centered around 0. If it
 * succeeds it returns true, if it fails because N is not a perfect square it returns false.
 */
bool grid(const int N, value_type * const x, value_type * const y, const value_type xRange, const value_type yRange){
    // Compute the number of gridpoints in one row or column
    int M = (int)sqrt(N);

    // Test if M is an integer
    if (N != M*M){
        std::cout << "The number of gridpoints is not a perfect square, use a different N." << std::endl;
        return false;
    }

    // Compute parameters
    const value_type dx = xRange/(M-1);
    const value_type dy = yRange/(M-1);
    const value_type xLim = xRange/2;
    const value_type yLim = yRange/2;

    // Loop over y coordinates
    for (int i = 0; i < M; ++i) {
        // Loop over x coordinates
        for (int j = 0; j < M; ++j) {
            x[i*M + j] = j*dx - xLim;
            y[i*M + j] = -i*dy + yLim;
        }
    }

    return true;
}


/*
 * This function creates a Lamb-Oseen vortex.
 */
void lambOseen(const int N, value_type * const x, value_type * const y, value_type * const q, const value_type visc, const value_type circ, value_type t, const value_type xCenter, const value_type yCenter){
    // Assuming the domain is centered around (0,0)
#pragma omp parallel for
    for (int i = 0; i < N; ++i) {
        // Compute the radius of the particle from the center of the vortex
        value_type r = sqrt((x[i]-xCenter)*(x[i]-xCenter) + (y[i]-yCenter)*(y[i]-yCenter));
        value_type sigma = sqrt(4*visc*(t+0.25/visc));  // If t = 0 then coreRadius sigma = 1
        q[i] = circ/(M_PI * sigma * sigma) * exp(-(r*r)/(sigma*sigma));
    }
}


/*
 * A function to smoothly cut-off some value array after a certain radius, using a cosine function
 */
void smoothCutOff(const int N, value_type * const q_in, value_type * const x_in, value_type * const y_in, const value_type radius_start, const value_type radius_end){
    // Loop over particles
    value_type radius;
    value_type radius_range = (radius_end - radius_start);
    value_type var;
#pragma omp parallel for
    for (int i = 0; i < N; ++i) {
        // Compute radius
        radius = sqrt(x_in[i]*x_in[i] + y_in[i]*y_in[i]);

        if(radius > radius_end){
            q_in[i] = 0;
        } else
        if(radius > radius_start){
            // Map [radius_start, radius_end] to [0, PI]
            var = (radius - radius_start)/radius_range * MPI;

            // Scale output with a cosine function of var
            q_in[i] = q_in[i] * (cos(var)+1)/2;

        }
        // If radius is smaller than radius_start do nothing.
    }
    return;
}


/*
 * This function copies the content of an array x into an Eigen MatrixXd.
 * Ordering:
 * Array - lexicographical (left to right then one row lower etc)
 * Matrix - filled column wise (so transpose of what data locations look like)
 */
void arrayToMatrix(value_type * x, MatrixXd & M, bool invertOrder){
    // Loop over rows and columns
    const int r = M.rows();
    const int c = M.cols();
    for (int i = 0; i < r; ++i) {
        for (int j = 0; j < c; ++j) {
            if(invertOrder){
                M(j,i) = x[i*c + j]; 
            }else{
                M(i,j) = x[i*c + j];
            }
        }
    }
    return;
}


/*
 * This function copies the content of an Eigen MatrixXd into an array.
 * Ordering:
 * Array - lexicographical
 * Matrix - filled column wise (so transpose of what data locations look like)
 */
void matrixToArray(value_type * x, MatrixXd & M, bool invertOrder){
    // Loop over rows and columns
    const int r = M.rows();
    const int c = M.cols();
    for (int i = 0; i < c; ++i) {
        for (int j = 0; j < r; ++j) {
            if(invertOrder){
                x[i*c + j] = M(j,i);
            } else {
                x[i * c + j] = M(i, j);
            }
        }
    }
    return;
}

/*
 * This function computes the analytical solution and writes the result to files.
 */
void analyticalSolution(){
    // PARAMETERS
    // time
    const value_type t_end = 0.1;
    const value_type dt = 0.0001;
    const int iter = t_end/dt;

    // vortex
    const value_type visc = 0.1;
    const value_type circ = 10;

    // spatial
    const int N = 40000;
    const int n =sqrt(N);
    const value_type dx = 0.01;
    const value_type range=(n-1)*dx;

    // Prepare string for filename
    std::string filenameX;
    std::string filenameY;
    std::string filenameQ;

    // Create a grid
    value_type * const x = new value_type[N];
    value_type * const y = new value_type[N];
    value_type * const q = new value_type[N];

    grid(N, x, y, range, range);


    // Loop over the timesteps and compute the vorticity field at each timestep. Write the resulting field to a file.
    for (int i = 0; i < iter; ++i) {
        lambOseen(N, x, y, q, visc, circ, i*dt, 0, 0);

        // Assign a string for filename and write to file
        filenameX = std::to_string(i) + "_X.txt";
        write_to_file(filenameX.c_str(), N, x);
        filenameY = std::to_string(i) + "_Y.txt";
        write_to_file(filenameY.c_str(), N, y);
        filenameQ = std::to_string(i) + "_Q.txt";
        write_to_file(filenameQ.c_str(), N, q);

    }

}
