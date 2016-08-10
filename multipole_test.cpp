
#include <iostream>
#include <chrono>
#include <cstring>
#include <cassert>
#include <cmath>
#include <cstdio>
#include "multipole.h"

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

int main() 
{
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

	// compute potential
	potential(theta_dist, xsrc, ysrc, sources, NSRC, xdst, ydst, NDST, xtargets);
	
	std::cout << "after POTENTIAL " << std::endl;
	
	// reference solution (from exercise solution)
	const int OFFSET = 0;
	const int JUMP = 1;

			//#pragma omp parallel for
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

	std::vector<double> a, b, c, d;

	for(int i = OFFSET; i < NDST; i += JUMP)
	{
		a.push_back(xref[i]);
		b.push_back(xtargets[i]);
		c.push_back(yref[i]);
		d.push_back(ytargets[i]);
	}

	check(&a[0], &b[0], a.size());
	
	// the errors obtained in the course solution are
	// l-infinity errors: 5.504e-09 (absolute) 1.713e-06 (relative)
    // l-1 errors: 2.597e-07 (absolute) 4.258e-05 (relative)

	
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

	
	return 0; 
	
}
