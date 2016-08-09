//
// Created by Isabelle Tan on 27-07-16.
//

#ifndef VORTEX_QUADTREE_H
#define VORTEX_QUADTREE_H

#include "simulation.h"			// for simulation parameters
#include "morton.h"

/*
 * The datastructure for nodes in the quad-Tree
 */
struct Node
{
    int level, morton_id;       // Nr of subdivisions in the tree for this node and morton index of square of this node.
    int child_id;               // Array index of the first of the four child nodes.
    int part_start, part_end;   // Array indices of first and last particle inside this node.
    value_type mass, xcom, ycom;// The mass of this node (sum of particle mass) and x,y position of center of mass.
    value_type r; 				// the radius of the smallest circle bounding the node with center at the center of mass of the node
    value_type* rxps; 			// real      part of p2e expansion
    value_type* ixps; 			// imaginary part of p2e expansion
    
};


/*
 * Build a tree
 */
void build(const value_type* const x, const value_type* const y, const value_type* const mass, const int N, const int k, value_type* xsorted, value_type* ysorted, value_type* mass_sorted, Node* tree, int depth);

/*
 * A function that splits a parent node recursively as long as there are more than k particles inside
 */
void split(Node* parent, Node* tree, int depth, int* index, value_type* xsorted, value_type* ysorted, value_type* mass_sorted, int k, int *newNodeIndex);

/*
 * Subdivide a parent node
 */
void assignParticles(Node* parent, Node* children, int depth, int* index);

/*
 * Compute the center of mass and the total mass in 4 children nodes.
 */
void centerOfMass(Node* children, value_type* xsorted, value_type* ysorted, value_type* mass_sorted);

/*
 * A function that prints all the attributes of a Node
 */
void printNode(Node node);


#endif //VORTEX_QUADTREE_H
