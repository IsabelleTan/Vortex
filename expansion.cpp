//
// Created by Isabelle Tan on 27-07-16.
//

#include "expansion.h"
#include "complexAVX.h"
#include <iostream>
#include <complex.h>
#include <limits>

/*
 * A kernel to compute the particle to expansion computations (simple version, without using AVX intrinsics)
 * This function computes the alpha_k coefficients and stores their real and imaginary
 * parts in the arrays rexpansion and iexpansion.
 */
void p2e_simple(const value_type* const x, const value_type* const y, const value_type* const q, const int nsources, const int order, const value_type xcom, const value_type ycom, value_type* const rexpansion, value_type* const iexpansion){
    // INPUT
    // x			array containing the x-coordinates of the source-particles 
    // y			.................... y-coordinates .......................
    // q			.................... the masses    .......................
    // nsources		number of sources 
    // order		order of the expansion
    // xcom			center of mass (x-coordinate)
    // ycom			.............. (y-coordinate)
    // OUTPUT
    // rexpansion	array of real 	   parts of the expansion 
    // iexpansion	........ imaginary ......................
    
    // Make arrays to store the radii and the values for different orders of r that are used in the multiplication
    std::complex<value_type> *r = new std::complex<value_type>[nsources];
    std::complex<value_type> *z = new std::complex<value_type>[nsources];

    // Fill the r array with initial values and fill the z array with the radii
    for (int i = 0; i < nsources; ++i) {
        r[i] = -q[i];
        z[i] = std::complex<value_type>(x[i], y[i]) - std::complex<value_type>(xcom, ycom);
    }

    // Compute the alpha's
    for (int j = 1; j <= order; ++j) {
        rexpansion[j - 1] = 0;
        iexpansion[j - 1] = 0;
        for (int k = 0; k < nsources; ++k) {
            // Multiply with z to obtain the next power of q*z^j
            r[k] *= z[k];

            // Add real and imaginary part to the arrays containing the expansion alpha's
            rexpansion[j - 1] += r[k].real();
            iexpansion[j - 1] += r[k].imag();
        }

        // Divide the whole sum by the order j
        rexpansion[j - 1] /= j;
        iexpansion[j - 1] /= j;
    }

    delete[] r;
	delete[] z;
    
};

/*
 * A kernel to compute the particle to expansion computations (using AVX intrinsics)
 * This function computes the alpha_k coefficients and stores their real and imaginary
 * parts in the arrays rexpansion and iexpansion.
 */
void p2e(const value_type* const x, const value_type* const y, const value_type* const q, const int nsources, const int order, const value_type xcom, const value_type ycom, value_type* const rexpansion, value_type* const iexpansion){
    // INPUT
    // x			array containing the x-coordinates of the source-particles 
    // y			.................... y-coordinates .......................
    // q			.................... the masses    .......................
    // nsources		number of sources 
    // order		order of the expansion
    // xcom			center of mass (x-coordinate)
    // ycom			.............. (y-coordinate)
    // OUTPUT
    // rexpansion	array of real 	   parts of the expansion 
    // iexpansion	........ imaginary ......................
    
    // Compute the vectorization parameters
    int SIMD_width = sizeof(__m256d) / sizeof(value_type);
    int SIMD_blocks = nsources / SIMD_width;

    // Make arrays to store the radii and the values for different orders of r so they can be reused in the computation
    value_type *r = new value_type[2 * nsources];
    value_type *z = new value_type[2 * nsources];

    // Fill the r array with initial values and fill the z array with the radii
    for (int i = 0; i < nsources; ++i) {
        r[i] = -q[i];
        r[nsources + i] = 0;
        z[i] = x[i] - xcom;
        z[nsources + i] = y[i] - ycom;
    }

    // Compute the alpha's using SIMD blocks
    for (int j = 1; j <= order; ++j) {
        rexpansion[j - 1] = 0;
        iexpansion[j - 1] = 0;
        for (int k = 0; k < SIMD_blocks; ++k) {
            // Load the values into registers
            __m256d x1 = _mm256_loadu_pd(r + k * SIMD_width);
            __m256d y1 = _mm256_loadu_pd(r + nsources + k * SIMD_width);
            __m256d x2 = _mm256_loadu_pd(z + k * SIMD_width);
            __m256d y2 = _mm256_loadu_pd(z + nsources + k * SIMD_width);
            __m256d res_r = _mm256_set1_pd(0);
            __m256d res_i = _mm256_set1_pd(0);


            // Multiply with z to obtain the next power of q*z^j
            multiply(x1, y1, x2, y2, &res_r, &res_i);

            // Store the new value in r
            _mm256_storeu_pd(r + k * SIMD_width, res_r);
            _mm256_storeu_pd(r + nsources + k * SIMD_width, res_i);

            // Add the real and imaginary part to the arrays containing the expansion alpha's
            rexpansion[j - 1] += HsumAvxDbl(res_r);
            iexpansion[j - 1] += HsumAvxDbl(res_i);
        }

        // Compute the left-over elements
        for (int i = SIMD_width * SIMD_blocks; i < nsources; ++i) {
            // Update r with the new power of q*z^j
            value_type r_r = r[i];
            value_type r_i = r[i + nsources];
            r[i] = r_r * z[i] - r_i * z[nsources + i];
            r[i + nsources] = r_r * z[nsources + i] + z[i] * r_i;

            // Add the new power of q*z^j to the real and imaginary parts of the expansion
            rexpansion[j - 1] += r[i];
            iexpansion[j - 1] += r[i + nsources];
        }

        // Divide the whole sum by the order j
        rexpansion[j - 1] /= j;
        iexpansion[j - 1] /= j;
    }

    delete[] r;
	delete[] z;

}

/*
 * This function computes the potential (streamfunction) value at the target particle located
 * at xtarget, ytarget.
 */
value_type e2p(const double xtarget, const double ytarget, const double q, const int order, const double* const rxps, const double* const ixps){
    // INPUT 
    // xtarget		target location (x-coordinate) 
    // ytarget		............... (y-coordinate)	
    // q			sum of masses of all source-particles
    // order		order of expansion
    // rxps			array containging the real 	    parts of expansion 
    // ixps			..................... imaginary ..................
    // OUTPUT
    // return value: 	potential at target location
    
    const value_type eps = std::numeric_limits<value_type>::epsilon();
    
    // Make an array to store the different orders of z_target
    std::complex<value_type> orders(1,0);
    std::complex<value_type> z_target(xtarget, ytarget);
    std::complex<value_type> streamFunction = q * log(z_target + eps);

    for (int i = 0; i < order; ++i) {
        orders *= (z_target + eps);
        streamFunction += std::complex<value_type>(rxps[i], ixps[i]) / orders;
    }

	return streamFunction.real();

};

/*
 * The particle to particle kernel. This function computes the interactions between nsources source particles and one target particle.
 */
value_type p2p(const value_type* const xsources, const value_type* const ysources, const value_type* const q, const int nsources, const value_type xtarget, const value_type ytarget){
	// INPUT
	// xsources		array containing location of source-particles (x-coordinate)
	// ysources		............................................. (y-coordinate)
	// q			................ masses   ...................
	// nsources		number of sources 
	// xtarget		target location (x-coordinate) 
	// ytarget		............... (y-coordinate)
	// OUTPUT
	// return value: 	potential at target location
      
    const value_type eps = std::numeric_limits<value_type>::epsilon();
    value_type streamfunction = 0;
    for (int i = 0; i < nsources; ++i) {
        streamfunction += q[i] * log( sqrt( (xsources[i] - xtarget)*(xsources[i] - xtarget) + (ysources[i] - ytarget)*(ysources[i] - ytarget) ) + eps);
    }

	return streamfunction;  
      
};

