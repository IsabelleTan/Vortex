
#include <iostream>
#include <chrono>
#include <cstring>
#include <cassert>
#include <cmath>
#include "multipole.h"
#include "datapoints.h"

double tol = 1e-8;
//double tol = 1e-2;

// from solution to exercises
void check(const double * ref, const double * res, const int N)
{
	double linf = 0, l1 = 0, linf_rel = 0, l1_rel = 0;

	for(int i = 0; i < N; ++i)
	{
		assert(!std::isnan(ref[i]));
		assert(!std::isnan(res[i]));

		const double err = ref[i] - res[i];
		const double maxval = std::max(fabs(res[i]), fabs(ref[i]));
		const double relerr = err/std::max(1e-6, maxval);

		if (fabs(relerr) >= tol && fabs(err) >= tol)
			printf("%d: %e ref: %e -> %e %e\n", i, res[i], ref[i], err, relerr);

//		assert(fabs(relerr) < tol || fabs(err) < tol);
//		assert(fabs(err) < tol);
//		assert(fabs(relerr) < tol);

		l1 += fabs(err);
		l1_rel += fabs(relerr);

		linf = std::max(linf, fabs(err));
		linf_rel = std::max(linf_rel, fabs(relerr));
	}

	printf("l-infinity errors: %.03e (absolute) %.03e (relative)\n", linf, linf_rel);
	printf("       l-1 errors: %.03e (absolute) %.03e (relative)\n", l1, l1_rel);
}

// reference solution from exercizes
void exercise_sol(double * xref, int OFFSET, int JUMP, double * xdst, double * ydst, int NDST, double * xsrc, double * ysrc, double * sources, int NSRC, double eps)
{
	for(int i = OFFSET; i < NDST; i += JUMP)
	{
		const double xd = xdst[i];
		const double yd = ydst[i];

		double s = 0;

		for(int j = 0; j < NSRC; ++j)
		{
			const double xr = xd - xsrc[j];
			const double yr = yd - ysrc[j];
			const double r2 = xr * xr + yr * yr;
			const double f  = fabs(r2) > eps;
			s += 0.5 * f * log(r2 + eps) * sources[j];
		}
		xref[i] = s;
	}
}

/*
 * A function to test the potential() function.
 */
void potential_test(){
	// open file etc.
	char filename[256];
	strcpy(filename, "diegoBinaryN400");		// now filename contains "diegoBinary ...."
	FILE * f = fopen(filename, "r");
	assert(f && sizeof(double) == sizeof(double));


	// read variables from file, memory-align them etc.
	int NSRC;
	fread(&NSRC, sizeof(int), 1, f);

	double *xsrc, *ysrc, *sources;
	posix_memalign((void **)&xsrc, 32, sizeof(double) * NSRC);
	posix_memalign((void **)&ysrc, 32, sizeof(double) * NSRC);
	posix_memalign((void **)&sources, 32, sizeof(double) * NSRC);

	fread(xsrc, sizeof(double), NSRC, f);
	fread(ysrc, sizeof(double), NSRC, f);
	fread(sources, sizeof(double), NSRC, f);

	int NDST;
	fread(&NDST, sizeof(int), 1, f);

	double *xdst, *ydst, *xref, *yref;
	posix_memalign((void **)&xdst, 32, sizeof(double) * NDST);
	posix_memalign((void **)&ydst, 32, sizeof(double) * NDST);
	posix_memalign((void **)&xref, 32, sizeof(double) * NDST);
	posix_memalign((void **)&yref, 32, sizeof(double) * NDST);

	fread(xdst, sizeof(double), NDST, f);
	fread(ydst, sizeof(double), NDST, f);
	fread(xref, sizeof(double), NDST, f);

	const double eps = std::numeric_limits<double>::epsilon() * 10;

	double *xtargets, *ytargets;
	posix_memalign((void **)&xtargets, 32, sizeof(double) * NDST);
	posix_memalign((void **)&ytargets, 32, sizeof(double) * NDST);

	printf("Testing %s with %d sources and %d targets (theta %.3e)...\n", "POTENTIAL", NSRC, NDST, theta_dist);

	// compute potential with nlogn method
	potential(theta_dist, xsrc, ysrc, sources, NSRC, xdst, ydst, NDST, xtargets);

	const int OFFSET = 0;
	const int JUMP = 1;

	// reference solution (choose the one from exercises, or our p2p implementation
	// potential_p2p(theta_dist, xsrc, ysrc, sources, NSRC, xdst, ydst, NDST, xref);
	exercise_sol(xref, OFFSET, JUMP, xdst, ydst, NDST, xsrc, ysrc, sources, NSRC, eps);

	std::vector<double> a, b, c, d;

	for(int i = OFFSET; i < NDST; i += JUMP)
	{
		a.push_back(xref[i]);	  // potential computed above
		b.push_back(xtargets[i]); // potential computed by our potential()
		c.push_back(yref[i]);
		d.push_back(ytargets[i]);
	}

	check(&a[0], &b[0], a.size());

	// free the memory
	free(xdst);
	free(ydst);

	free(xtargets);
	free(ytargets);

	free(xref);
	free(yref);

	free(xsrc);
	free(ysrc);
	free(sources);

	return;
}

/*
 * A function to time the potential function.
 */
value_type potential_time(int N, value_type theta){
	// Prepare time variables
	auto start = std::chrono::high_resolution_clock::now();
	auto end   = std::chrono::high_resolution_clock::now();
	value_type potentialTime;

	// Prepare arrays for potential
	value_type * xsrc;
	value_type * ysrc;
	value_type * qsrc;
	value_type * xdst;
	value_type * ydst;
	value_type * potdst;

	// Align arrays
	posix_memalign((void **)&xsrc, 32, sizeof(double) * N);
	posix_memalign((void **)&ysrc, 32, sizeof(double) * N);
	posix_memalign((void **)&qsrc, 32, sizeof(double) * N);
	posix_memalign((void **)&xdst, 32, sizeof(double) * N);
	posix_memalign((void **)&ydst, 32, sizeof(double) * N);
	posix_memalign((void **)&potdst, 32, sizeof(double) * N);

	// Initialize arrays
	load_data_random(N, &xsrc, &ysrc);
	load_data_random(N, &qsrc);
	load_data_random(N, &xdst, &ydst); // (Note that here xdst and ydst are not on a grid)

	// Run simulation and time
	start = std::chrono::high_resolution_clock::now();
	potential(theta, xsrc, ysrc, qsrc, N, xdst, ydst,N, potdst);
	end = std::chrono::high_resolution_clock::now();
	potentialTime = static_cast<std::chrono::duration<double>>(end-start).count();

	delete[] xsrc;
	delete[] ysrc;
	delete[] qsrc;
	delete[] xdst;
	delete[] ydst;
	delete[] potdst;

	return potentialTime;
}

/*
 * A function to time something N times and compute the average and variance
 */
void time(int N){
	value_type *times = new value_type[N];

	// Print parameters
	std::cout << "Timing " << N << " simulations of potential(), with" << std::endl;
	std::cout << "nParticles = " << nParticles << "\ntheta_dist = " << theta_dist << "\nexp_order = " << exp_order << "\n" << std::endl;

	// Track the time N times
	for (int i = 0; i < N; ++i) {
		std::cout << "Timing " << i << " of " << N << std::endl;
		times[i] = potential_time(nParticles, theta_dist);
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
	std::cout << "\nTimed " << N << " times." << std::endl;
	std::cout << "Mean = " << mean << std::endl;
	std::cout << "Variance = " << var << " is " << var/mean * 100 << "% of the mean" << std::endl;

	delete[] times;
}

int main()
{
	time(100);
	//potential_test();
	
	return 0; 
	
}

/* RESULTS
 *
 * Comparing our p2p to exercize solution:
 * l-infinity errors: 1.153e-01 (absolute) 1.884e+00 (relative)
 * l-1 errors: 2.120e+00 (absolute) 1.457e+02 (relative)
 *
 * Comparing our potential to exercize solution
 * l-infinity errors: 1.153e-01 (absolute) 1.884e+00 (relative)
 * l-1 errors: 2.120e+00 (absolute) 1.457e+02 (relative)
 *
 * Comparing our potential to our p2p
 * l-infinity errors: 1.230e-08 (absolute) 5.733e-04 (relative)
 * l-1 errors: 6.798e-07 (absolute) 8.700e-04 (relative)
 *
 */


