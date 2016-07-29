
#include "morton.h"

// Put some random values in the arrays
void initialize(const int N, value_type* x, value_type* y, value_type* mass){
    for (int i = 0; i < N; ++i) {
        x[i] = drand48();
        y[i] = drand48();
        mass[i] = drand48();
    }
}

void extent(const int N, const double* const x, const double* const y, double& xmin, double& ymin, double& ext) 								// set 2, question 1, a)
{
	// N	# of data points
	// x	array of x-coordinates
	// y	array of y-coordinates
	// xmin min of all x-values
	// ymin	min of all y-values
	// ext 	the "extent", i.e. the size of smallest square including all data points
	
	// compute the min and max along the x- and y-axis
	// launch on 4 threads
	std::thread t1([&N, &x, &xmin]() mutable{			
		xmin = *(std::min_element(&(x[0]), &(x[N]))); 					// find the min along the x-axis and assign it to xmin
	}); 

	std::thread t2([&N, &y, &ymin]() mutable{							
		ymin = *(std::min_element(&(y[0]), &(y[N])));					// find the min along the y-axis and assign it to ymin
	}); 

	double xmax;
	std::thread t3([&N, &x, &xmax]() mutable{							
		xmax = *(std::max_element(&(x[0]), &(x[N])));					// find the max along the x-axis and assign in to xmax
	}); 

	double ymax = *(std::max_element(&(y[0]), &(y[N])));				// find the max along the y-axis and assign in to xmin
																		// this runs on the master thread
	t1.join();															// wait for the slave-threads to complete their task and join them to the master
	t2.join();
	t3.join();
	
	// compute the extent and shift the domain boundaries for numerical reasons
	const double eps = 10000 * std::numeric_limits<double>::epsilon();	
	ext = std::max(xmax - xmin, ymax - ymin) + eps;						
	xmin -= eps * ext; 
	ymin -= eps * ext; 
} 
	
void morton(const int N, const double* const x, const double* const y, const double xmin, const double ymin, const double ext, int* index) 		// set 2, question 1, b)
{
	// N		# of data points
	// x		array of x-coordinates 										(unsorted)
	// y		array of y-coordinates 										(unsorted)
	// index	array of Morton indices corresponding to the data points	(unsorted)
	
	// We don't care about recursive subdivision or anything like that
	// rather just assign particles their morton index in a 16-level tree [from 0 to 15]	
	// highest level = level 15 -> 2*15 = 30 bits 
	//               is defined in #define level_max at top of document
	
	/// problem: the particles are not sorted like a downwards Z, but rather like: 
	///							 _________
	///							|	 |    |
	///							| 10 | 11 | 
	///							|____|____|
	///							|	 |    |
	///							| 00 | 01 |
	///							|____|____| 
	///
		
	#pragma omp parallel for
    for(int i = 0; i < N; ++i)
    {
		int xid = floor((x[i] - xmin) / ext * (1 << depthtree));	// (x[i] - xmin) / ext // normalize the coord. pos. w/ resp. to the square length
		int yid = floor((y[i] - ymin) / ext * (1 << depthtree));	// * (1 << level_max)  // 0..01000000000000000 = 2^15

		xid = (xid | (xid << 8)) & 0x00FF00FF; // cf https://graphics.stanford.edu/~seander/bithacks.html#InterleaveBMN
		xid = (xid | (xid << 4)) & 0x0F0F0F0F;
		xid = (xid | (xid << 2)) & 0x33333333;
		xid = (xid | (xid << 1)) & 0x55555555;

		yid = (yid | (yid << 8)) & 0x00FF00FF;
		yid = (yid | (yid << 4)) & 0x0F0F0F0F;
		yid = (yid | (yid << 2)) & 0x33333333;
		yid = (yid | (yid << 1)) & 0x55555555;

		index[i] = xid | (yid << 1);			// assign morton index
    }
}

void sort(const int N, int* index, int* keys) 																									// set 2, question 1, c)
{
	// INPUT
	// N		# of data points
	// index	array of Morton indices corresponding to the data points	(unsorted)
	// keys		= {0, 1, 2, ..., N-1}
	// OUTPUT 
	// index	is now sorted in ascending order
	// keys		contains permutation that sorts "index".
	
	/// memory bandwidth is the limiting ressource
	
	/// can I avoid the "keys" and "temp" array 
	/// by just putting the x and y in a 3-uple with the indices and sorting that? 
	/// or is it actually more memory intensive cause you move the x and y elements around as well? 
	/// I think so ... 
		
	std::pair<int, int> *temp = NULL;						// create an array of pairs that will contain (1) index of a particle (2) key corresponding to its position in the x-array and y-array
	posix_memalign((void **)&temp, 32, sizeof(*temp)*N );	// memalign it 

	#pragma omp parallel for	
	for(size_t i=0; i<N; ++i){								// fill it up with values of indices and keys
		temp[i].first = index[i]; 
		temp[i].second = keys[i]; 
	}

	__gnu_parallel::sort(temp, temp+N);						// sort by ascending indices. Alongside this, the keys get automatically sorted as well 

	#pragma omp parallel for	
	for(size_t i=0; i<N; ++i){
		index[i] = temp[i].first;							/// do I actually need to reorder "index"? It is not reused in "void reorder"... is it reused somewhere else? If not, I can probably delete this line.
		keys[i] = temp[i].second;
	}

	free(temp);												// free the memory assigned by posix_memalign
}

void reorder(const int N, const int* const keys, const double* const x, const double* const y, const double* const q, double* xsorted, double* ysorted, double* qsorted) 				// set 2, question 1, d)
{
	// INPUT
	// N		# of data points
	// x		array of x-coordinates 									(unsorted)
	// y		array of y-coordinates 									(unsorted)
	// keys		contains permutation that sorts the particles according to their "index"
	// OUTPUT 
	// xsorted  array of x-coordinates sorted wrt keys
	// ysorted  ........ y-coordinates ...............
	// qsorted  ........ masses        ...............
	// 
	
	/// why do we refill a different array (xsorted, ysorted) instead of just sorting the arrays x and y themselves directly? 
	/// because then we might re-write over some values etc. It's safer and quicker to just write a new array

	#pragma omp parallel for
	for(size_t i=0; i<N; i++){
		const int key(keys[i]);		// pre-fetch value so I don't have to do a random access in the "keys"-array twice.
		xsorted[i] = x[key];
		ysorted[i] = y[key];
		qsorted[i] = q[key];
	}
}


