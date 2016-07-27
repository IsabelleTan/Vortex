//
// Created by Isabelle Tan on 27-07-16.
//

#ifndef VORTEX_QUADTREE_H
#define VORTEX_QUADTREE_H

// Set the value_type
typedef double value_type;

/*
 * The datastructure for nodes in the quad-Tree
 */
struct Node
{
    int level, morton_id;       // Nr of subdivisions in the tree for this node and morton index of square of this node.
    int child_id;               // Array index of the first of the four child nodes.
    int part_start, part_end;   // Array indices of first and last particle inside this node.
    value_type mass, xcom, ycom;     // The mass of this node (sum of particle mass) and x,y position of center of mass.
};

/*
 * The function that builds the quadtree
 */
void build(const value_type* const x, const value_type* const y, const value_type* const mass, const int N, const int k, value_type* xsorted, value_type* ysorted, value_type* mass_sorted, Node* tree, int depth);

/*
 * A function to print the content of a node.
 */
void printNode(Node node);


#endif //VORTEX_QUADTREE_H
