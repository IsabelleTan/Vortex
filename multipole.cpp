
#include "multipole.h"							

/*
 * This function computes the potential (i.e. the stream function psi) at all target locations (on a grid),
 * using the vorticity ("mass") of all source particles, which could be located anywhere in domain. 
 * We're using a recursive approach and starting with the root node
 */
//! parallel with respect to targets 
//! but load imbalance: traversal of tree is different for each target 
void potential(double theta, 
			  double *xsrc, double *ysrc, double *qsrc, int nsrc,
			  double *xdst, double *ydst, int ndst, double *potdst)
{
	// INPUT
	// theta	parameter used to define which nodes are considered as far neighbors from a certain 2D-location
	// xsrc		array containing the x-coordinates of the source-particles (all)
	// ysrc		.................... y-coordinates .......................
	// qsrc		.................... "masses" (i.e. vorticities) of the source-particles ("source-values")
	// nsrc		number of source particles
	// xdst		array containing the x-coordinates of the target-particles ("destination-particles") (all)
	// ydst		.................... y-coordinates .......................
	// ndst		number of target-particles
	// OUPUT
	// potdst	values of the potential at destination
	
	// create all the variables I will need
	const double thetasquared = theta*theta;
	double *xsorted, *ysorted, *qsorted; 
	Node *nodes;
	double *result, xt, yt;
	int maxnodes = (nsrc + kleaf - 1) / kleaf * 60;

	// do the memory-alignments
	posix_memalign((void **)&xsorted, 32, sizeof(double) * nsrc);
	posix_memalign((void **)&ysorted, 32, sizeof(double) * nsrc);
	posix_memalign((void **)&qsorted, 32, sizeof(double) * nsrc);
	posix_memalign((void **)&nodes,   32, sizeof(Node) * maxnodes);
	
	// build the quadtree
	build(xsrc, ysrc, qsrc, nsrc, kleaf, xsorted, ysorted, qsorted, nodes, depthtree);
		
	// evaluate the potential field at the location of each target point:
	for(size_t i=0; i<ndst; i++){
		
		result = 0;
		
		//! naive: start at level 0:
		//! evaluate(nodes, 0, xsorted, ysorted, qsorted, thetasquared, result, xdst[i], ydst[i]);
		//! but I know there are no leaves and most probably no squares far away enough in the two first levels (level 0 and level 1) 
		//! so I might as well spare some computation and start on level 2 
		//! there are 16 squares on level 2:
		
		evaluate(nodes,  0, xsorted, ysorted, qsorted, thetasquared, result, xdst[i], ydst[i]); /* square # 00 00 */
		evaluate(nodes,  1, xsorted, ysorted, qsorted, thetasquared, result, xdst[i], ydst[i]); /* square # 00 01 */
		evaluate(nodes,  2, xsorted, ysorted, qsorted, thetasquared, result, xdst[i], ydst[i]); /* square # 00 10 */
		evaluate(nodes,  3, xsorted, ysorted, qsorted, thetasquared, result, xdst[i], ydst[i]); /* square # 00 11 */
		
		evaluate(nodes,  4, xsorted, ysorted, qsorted, thetasquared, result, xdst[i], ydst[i]); /* square # 01 00 */
		evaluate(nodes,  5, xsorted, ysorted, qsorted, thetasquared, result, xdst[i], ydst[i]); /* square # 01 01 */
		evaluate(nodes,  6, xsorted, ysorted, qsorted, thetasquared, result, xdst[i], ydst[i]); /* square # 01 10 */
		evaluate(nodes,  7, xsorted, ysorted, qsorted, thetasquared, result, xdst[i], ydst[i]); /* square # 01 11 */
		
		evaluate(nodes,  8, xsorted, ysorted, qsorted, thetasquared, result, xdst[i], ydst[i]); /* square # 10 00 */
		evaluate(nodes,  9, xsorted, ysorted, qsorted, thetasquared, result, xdst[i], ydst[i]); /* square # 10 01 */
		evaluate(nodes, 10, xsorted, ysorted, qsorted, thetasquared, result, xdst[i], ydst[i]); /* square # 10 10 */
		evaluate(nodes, 11, xsorted, ysorted, qsorted, thetasquared, result, xdst[i], ydst[i]); /* square # 10 11 */
		
		evaluate(nodes, 12, xsorted, ysorted, qsorted, thetasquared, result, xdst[i], ydst[i]); /* square # 11 00 */
		evaluate(nodes, 13, xsorted, ysorted, qsorted, thetasquared, result, xdst[i], ydst[i]); /* square # 11 01 */
		evaluate(nodes, 14, xsorted, ysorted, qsorted, thetasquared, result, xdst[i], ydst[i]); /* square # 11 10 */
		evaluate(nodes, 15, xsorted, ysorted, qsorted, thetasquared, result, xdst[i], ydst[i]); /* square # 11 11 */
		
		potdst[i] = *result;  	// assign result 
	}
	
	// free the memory allocated
	free(xsorted);
	free(ysorted);
	free(qsorted);
	free(nodes);
	
}


/*
 * Evaluate the potential (i.e. the stream function psi) at a particular target location (xt, yt), 
 * given an array of nodes and their expansions 
 */
void evaluate(const Node* nodes, const int node_id, /*const double* expansions // putting it as attribute of class Node instead, */ const double *xdata, const double *ydata, const double *mdata,
		const double thetasquared, double * const result, const double xt, const double yt)
{
	// INPUT
	// nodes	array of nodes for the tree of data we're working with
	// node_id	index indicating the node I am considering now (used for recursive call) inside the array of nodes 
	// expansions	array of expansions for these nodes 
	// xdata	array of x-coordinates of the source-particles (all of them, not just the ones contained in the node. Needed for passing this array down in the recursion)
	// ydata	........ y-coordinates .......................
	// mdata	........ source-values ....................... (vorticities)
	// thetasquared	.
	// xt		x-coordinate of the location of the target
	// yt		y-coordinate .............................
	// OUTPUT
	// result	potential at target location (xt, yt)
	
	const Node* const node = nodes + node_id;														// create a pointer to the node we're considering
	double distsquared = (node->xcom - xt)*(node->xcom - xt) + (node->ycom - yt)*(node->ycom -yt) ; // square of the distance between the target location and the center of mass of the node 
	
	// recursive algorithm for the evaluation of the potential 
	//! should I add a factor 4 here under as in the pseudo-code? (lecture 4, slide 24)
	if( (node->r * node->r) < (distsquared * thetasquared) ){			// if the target and the node are far away from each other
																		// return the expansion to particle expression as a potential
		*result = e2p(xt - node->xcom, yt - node->ycom, node->mass, exp_order, node->rxps, node->ixps);	
		return;
	} else if(node->child_id == -1){ 									// if this node is a leaf, we cannot descend deeper into the tree
																		// return the particle-to-particle expansion
		const int node_start = node->part_start;
		*result = p2p(&xdata[node_start], &ydata[node_start], &mdata[node_start], node->part_end - node->part_start, xt, yt);
		return;
	} else {															// in any other case, continue the descent of the tree
		//! should I set result to zero here? As stated in exercizes? But then the variable result will be set back to zero at each recursive call, seems weird
		evaluate(nodes, node->child_id  , xdata, ydata, mdata, thetasquared, result, xt, yt) ;
		evaluate(nodes, node->child_id+1, xdata, ydata, mdata, thetasquared, result, xt, yt) ;
		evaluate(nodes, node->child_id+2, xdata, ydata, mdata, thetasquared, result, xt, yt) ;
		evaluate(nodes, node->child_id+3, xdata, ydata, mdata, thetasquared, result, xt, yt) ;
	}	
}


