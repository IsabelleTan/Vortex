
#include <iostream>							
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
	
	std::cout << "in POTENTIAL " << std::endl;

	// create all the variables needed
	const double thetasquared = theta*theta;
	double *xsorted, *ysorted, *qsorted; 
	Node *nodes;
	double *result(new double(0));
	int maxnodes = (nsrc + kleaf - 1) / kleaf * 60;


	// do the memory-alignments
	posix_memalign((void **)&xsorted, 32, sizeof(double) * nsrc);
	posix_memalign((void **)&ysorted, 32, sizeof(double) * nsrc);
	posix_memalign((void **)&qsorted, 32, sizeof(double) * nsrc);
	posix_memalign((void **)&nodes,   32, sizeof(Node) * maxnodes);

	
	// build the quadtree
	build(xsrc, ysrc, qsrc, nsrc, kleaf, xsorted, ysorted, qsorted, nodes, depthtree);
	
	// evaluate the potential field at the location of each target point:
	//! 	#pragma omp parallel for schedule(static,1)
	for(size_t i=0; i<ndst; i++){

		// for one target point (xdst[i], ydst[i]), do the following
		*result = 0;

		//! the naive approach would be to start at level 0:
		//! evaluate(nodes, 0, xsorted, ysorted, qsorted, thetasquared, result, xdst[i], ydst[i]);
		//! but we know there are no leaves and most probably no squares far away enough in the two first levels (level 0=root and level 1) 
		//! so we might as well spare ourselves some computation and start on level 2 
		//! there are 16 squares on level 2:
		
		// nodes[0] 	root node	level 0 
		// nodes[1] 	child 0 	level 1
		// nodes[2] 	child 1	 
		// nodes[3] 	child 2	 
		// nodes[4] 	child 3	 
		// nodes[5] 	...			level 2 	 
		
		// get the indices of the 16 level-2 nodes 
		// 5 = nodes[1]->child_id			level 2, children of node # 1
		// 6
		// 7
		// 8
		// s2 = nodes[2]->child_id			level 2, children of node # 2
		// s2+1
		// s2+2
		// s2+3
		
		unsigned int s2 = nodes[2].child_id;
		unsigned int s3 = nodes[3].child_id;
		unsigned int s4 = nodes[4].child_id;
				
		evaluate(nodes,    5, xsorted, ysorted, qsorted, thetasquared, result, xdst[i], ydst[i]); /* square # 00 00 */ 	//std::cout << "done with level 2.0  " << std::endl;
		evaluate(nodes,    6, xsorted, ysorted, qsorted, thetasquared, result, xdst[i], ydst[i]); /* square # 00 01 */
		evaluate(nodes,    7, xsorted, ysorted, qsorted, thetasquared, result, xdst[i], ydst[i]); /* square # 00 10 */
		evaluate(nodes,    8, xsorted, ysorted, qsorted, thetasquared, result, xdst[i], ydst[i]); /* square # 00 11 */
		
		evaluate(nodes, s2  , xsorted, ysorted, qsorted, thetasquared, result, xdst[i], ydst[i]); /* square # 01 00 */	//std::cout << "hello 6 " << std::endl;
		evaluate(nodes, s2+1, xsorted, ysorted, qsorted, thetasquared, result, xdst[i], ydst[i]); /* square # 01 01 */
		evaluate(nodes, s2+2, xsorted, ysorted, qsorted, thetasquared, result, xdst[i], ydst[i]); /* square # 01 10 */ //std::cout << "hello 6 " << std::endl;
		evaluate(nodes, s2+3, xsorted, ysorted, qsorted, thetasquared, result, xdst[i], ydst[i]); /* square # 01 11 */
		
		evaluate(nodes, s3  , xsorted, ysorted, qsorted, thetasquared, result, xdst[i], ydst[i]); /* square # 01 00 */	//std::cout << "hello 7 " << std::endl;
		evaluate(nodes, s3+1, xsorted, ysorted, qsorted, thetasquared, result, xdst[i], ydst[i]); /* square # 01 01 */
		evaluate(nodes, s3+2, xsorted, ysorted, qsorted, thetasquared, result, xdst[i], ydst[i]); /* square # 01 10 */
		evaluate(nodes, s3+3, xsorted, ysorted, qsorted, thetasquared, result, xdst[i], ydst[i]); /* square # 01 11 */
		
		evaluate(nodes, s4  , xsorted, ysorted, qsorted, thetasquared, result, xdst[i], ydst[i]); /* square # 01 00 */	//std::cout << "hello 8 " << std::endl;
		evaluate(nodes, s4+1, xsorted, ysorted, qsorted, thetasquared, result, xdst[i], ydst[i]); /* square # 01 01 */
		evaluate(nodes, s4+2, xsorted, ysorted, qsorted, thetasquared, result, xdst[i], ydst[i]); /* square # 01 10 */
		evaluate(nodes, s4+3, xsorted, ysorted, qsorted, thetasquared, result, xdst[i], ydst[i]); /* square # 01 11 */
			//std::cout << "hello 18 " << std::endl;
		potdst[i] = *result;  	// assign result 
			//std::cout << "hello 19 " << std::endl;

	}
	
	// free the memory allocated
	free(xsorted);
	free(ysorted);
	free(qsorted);
	free(nodes);
	delete result;
	
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
	// node_id	index indicating the node we am considering now (used for recursive call) inside the array of nodes 
	// expansions	array of expansions for these nodes (not needed, we're actually implementing expansion as an attribute of node)
	// xdata	array of x-coordinates of the source-particles (all of them, not just the ones contained in the node. Needed for passing this array down in the recursion)
	// ydata	........ y-coordinates .......................
	// mdata	........ source-values ....................... (vorticities)
	// thetasquared	.
	// xt		x-coordinate of the location of the target
	// yt		y-coordinate .............................
	// OUTPUT
	// result	potential at target location (xt, yt)
	
//	std::cout << "in  EVALUATE " << std::endl;
	
	const Node* const node = nodes + node_id;														// create a pointer to the node we're considering
//	std::cout << "got node  " << std::endl;
//	printNode(*node);

	if(node->r > -0.5){ // if this node isn't empty, go on with the evaluation the whole evaluation and stop the descent
		
		if(node->r == 0){ // if the node contains only 1 particle 
			// compute the p2p expansion
			const int n_s = node->part_start;
			*result += p2p(xdata + n_s, ydata + n_s, mdata + n_s, 1, xt, yt);
			return ; 
		}
		
		// if the node contains strictly more than one particle, compute the distance between the target and the center of mass of the node 
		double distsquared = ( (node->xcom - xt)*(node->xcom - xt) + (node->ycom - yt)*(node->ycom -yt) ) ; // square of the distance between the target location and the center of mass of the node 

		// recursive algorithm for the evaluation of the potential 
		// pseudo-code, cf lecture 4, slide 24
		if( (node->r * node->r) < (distsquared * thetasquared) ){			// if the target and the node are sufficiently far away from each other
																			// return the expansion to particle expression as a potential
//			std::cout << "e2p expansion" << std::endl;
			*result += e2p(xt - node->xcom, yt - node->ycom, node->mass, exp_order, node->rxps, node->ixps);	
			return;
		} else if(node->child_id == -1){ 									// if they're not sufficiently far away from each other
																			// and if this node is a leaf, 
																			// then we cannot descend deeper into the tree
																			// return the particle-to-particle expansion
//			std::cout << "p2p expansion" << std::endl;
			const int n_s = node->part_start;
			*result += p2p(xdata + n_s, ydata + n_s, mdata + n_s, node->part_end - n_s, xt, yt);
			return;
		} else {															// in any other case, 
																			// i.e., they are not sufficiently far away from each other, and this node isn't a leaf
																			// continue the descent of the tree
//			std::cout << "descending into tree" << std::endl;
			evaluate(nodes, node->child_id  , xdata, ydata, mdata, thetasquared, result, xt, yt) ;
			evaluate(nodes, node->child_id+1, xdata, ydata, mdata, thetasquared, result, xt, yt) ;
			evaluate(nodes, node->child_id+2, xdata, ydata, mdata, thetasquared, result, xt, yt) ;
			evaluate(nodes, node->child_id+3, xdata, ydata, mdata, thetasquared, result, xt, yt) ;
		}
		
	}
	
	*result += 0;	

}
