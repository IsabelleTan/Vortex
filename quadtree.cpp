//
// Created by Isabelle Tan on 27-07-16.
//

#include <iostream>
#include "quadtree.h"

using namespace std;


/*
 * The function that builds the quadtree.
 */
void build(const value_type* const x, const value_type* const y, const value_type* const mass, const int N, const int k, value_type* xsorted, value_type* ysorted, value_type* mass_sorted, Node* tree, int depth){
    // TODO complete build
}

/*
 * A function that prints the attributes of an object of type Node
 */
void printNode(Node node){
    cout<< "\n\nNode at address: " << &node << endl;
    cout<< "level: " << node.level << endl;
    cout<< "morton_id: " << node.morton_id << endl;
    cout<< "child_id: " << node.child_id << endl;
    cout<< "part_start: " << node.part_start << endl;
    cout<< "part_end: " << node.part_end << endl;
    cout<< "mass: " << node.mass << endl;
    cout<< "xcom: " << node.xcom << endl;
    cout<< "ycom: " << node.ycom << endl;
}