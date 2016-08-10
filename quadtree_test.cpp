//
// Created by Isabelle Tan on 04-05-16.
//

#include <iostream>
#include "morton.h"
#include "quadtree.h"

// Set the value_type
typedef double value_type;
value_type epsilon = 0.000001;

// Test the initialization function
bool testInitialization(){
    int N = 100;

    // Set return bool to true
    bool result = true;

    // Some test particle coordinates
    value_type *x = new value_type[N];
    value_type *y = new value_type[N];
    value_type *mass = new value_type[N];

    // Put some values in the arrays
    initialize(N, x, y, mass);

    // Check whether the x and y values are within the range (0,1)
    for (int i = 0; i < N; ++i) {
        if(x[i] < 0 || x[i] > 1 || y[i] < 0 || y[i] > 1){
            result = false;
        }
    }

    // Release the allocated memory
    delete[] x;
    delete[] y;

    return result;
}

// Test the computation of the extent and location (xmin,ymin)
bool testComputeExtent(){
    int N = 2;
    // Set the resulting boolean to true
    bool result = true;

    // Initiate arrays
    value_type *x = new value_type[N];
    value_type *y = new value_type[N];
    x[0] = 0.1;
    x[1] = 0.9;
    y[0] = 0.05;
    y[1] = 0.85;

    // Initiate extend and location variables
    value_type xmin;
    value_type ymin;
    value_type ext;

    // Compute the extend and xmin, ymin
    extent(2, x, y, xmin, ymin, ext);

    // Check whether the values are correct
    if(xmin != x[0] || ymin != y[0] || ext - 0.8 > epsilon){
        result = false;
    }

    // Release the allocated memory
    delete[] x;
    delete[] y;

    return result;
}

// Test the function of center of mass
bool testCenterOfMass(){
    bool result = true;
    // Make lists of masses and xy coordinates
    int N = 7;

    value_type xsorted[N] = {0.2, 0.25, 0.4, 0.48, 0.7, 0.72, 0.8};
    value_type ysorted[N] = {0.3, 0.2, 0.2, 0.1, 0.7, 0.4, 0.4};
    value_type masssorted[N] = {0.4, 0.3, 0.3, 0.89, 0.41, 0.1, 0.66};

    // Initiate the children
    Node children[4];

    // Initialize the children nodes
    Node child_0 = Node { 1, // level
                          0, // morton index
                          -1, // child_id
                          0, // part_start
                          2, // part_end
                          1, // node mass
                          1, // x center of mass
                          1  // y center of mass
    };

    Node child_1 = Node { 1, // level
                          0, // morton index
                          -1, // child_id
                          3, // part_start
                          4, // part_end
                          1, // node mass
                          1, // x center of mass
                          1  // y center of mass
    };

    Node child_2 = Node { 1, // level
                          0, // morton index
                          -1, // child_id
                          5, // part_start
                          5, // part_end
                          1, // node mass
                          1, // x center of mass
                          1  // y center of mass
    };

    Node child_3 = Node { 1, // level
                          0, // morton index
                          -1, // child_id
                          6, // part_start
                          6, // part_end
                          1, // node mass
                          1, // x center of mass
                          1  // y center of mass
    };

    // Put the children nodes into the array
    children[0] = child_0;
    children[1] = child_1;
    children[2] = child_2;
    children[3] = child_3;

    centerOfMass(children, xsorted, ysorted, masssorted);

    // Test the values against known true values
    value_type xcom_true[4] = {0.275, 0.5493846, 0.72, 0.8};
    value_type ycom_true[4] = {0.24, 0.2892307, 0.4, 0.4};
    value_type mass_true[4] = {1, 1.3, 0.1, 0.66};

    for (int i = 0; i < 4; ++i) {
        if (children[i].mass - mass_true[i] > 0.000001) {
            std::cout << "Mass of node " << i << " computed incorrectly." << std::endl;
            result = false;
        }
        if (children[i].xcom - xcom_true[i] > 0.000001) {
            std::cout << "xcom of node " << i << " computed incorrectly." << std::endl;
            result = false;
        }
        if (children[i].ycom - ycom_true[i] > 0.000001) {
            std::cout << "ycom of node " << i << " computed incorrectly." << std::endl;
            result = false;
        }
    }

    if (result) {
        std::cout << "Test succeeded." << std::endl;
    }
    return result;
}

// Test the split function
bool testBuild(){
    // Test the stopping criterion (#p=7, k=8).
    bool result = true;

    int N = 7;
    int k = 8;

    value_type x[7] = {12,45,34,90,34,23,56};
    value_type y[7] = {16,82,72,2,45,89,52};
    value_type mass[7] = {1, 2, 3, 4, 5, 6, 7};

    value_type xsorted[7];
    value_type ysorted[7];
    value_type mass_sorted[7];

    // Allocate the tree array containing all the nodes
    int maxNodes = (int) std::min((float)8 * N / k, (float)pow(4, depthtree));
    Node* tree = new Node[maxNodes];

    build(x, y, mass, N, k, xsorted, ysorted, mass_sorted, tree, depthtree);

    // Test the resulting tree
    if (tree[0].child_id != -1) {
        std::cout << "Test of stopping criterion failed, the root node is not a leaf node." << std::endl;
        std::cout << "root.child_id = " << tree[0].child_id << std::endl;
        result = false;
    } else {
        std::cout << "Root node is an empty node => The test for the stopping criterion passed." << std::endl;
    }


    // Test building the tree for k=2
    N=14;
    k=2;

    // Create new particles
    value_type x2[14] = {12,45,34,90,34,23,56,3,76,54,32,56,45,90};
    value_type y2[14] = {16,82,72,2,45,89,52,54,12,12,12,43,2,89};
    value_type mass2[14] = {1, 0.2, 3, 12, 5, 6, 3, 13, 9, 2, 11, 12, 3, 14};

    // Allocate arrays for sorted particles
    value_type xsorted2[14];
    value_type ysorted2[14];
    value_type mass_sorted2[14];

    // Compute new maxNodes for treesize
    int maxNodes2 = (int) std::min((float)8 * N / k, (float)pow(4, depthtree));

    Node* tree2 = new Node[maxNodes2];

    std::cout << "Testing the build() function for a tree with N=14, k=2 and depth=16." << std::endl;
    build(x2, y2, mass2, N, k, xsorted2, ysorted2, mass_sorted2, tree2, depthtree);

    // Test the resulting tree
    for (int i = 0; i < maxNodes2; ++i) {
        printNode(tree2[i]);
    }

    if (result) {
        std::cout << "Check the created nodes to validate the function." << std::endl;
    }

    delete[] tree;
    delete[] tree2;
    return result;
}

// Test the computation of attribute "r"
bool testr(){
	
	bool result = true; 
	
	// Create particles and particle arrays 
	int N = 10;
    int k = 8;

    value_type x[10] = {0.1, 0.15, 0.4, 0.6, 0.7, 0.71, 0.8, 0.74, 0.38, 0.41};
    value_type y[10] = {0.6, 0.85, 0.65, 0.3, 0.6, 0.1, 0.4, 0.33, 0.42, 0.57};
    value_type mass[10] = {1, 2, 3, 4, 5, 6, 7, 2.3, 1.4, 5.7};

    value_type xsorted[10];
    value_type ysorted[10];
    value_type mass_sorted[10];

    // Allocate the tree array containing all the nodes
    int maxNodes = (int) std::min((float)8 * N / k, (float)pow(4, depthtree));
    Node* tree = new Node[maxNodes];

	// create a tree to be tested 
    build(x, y, mass, N, k, xsorted, ysorted, mass_sorted, tree, depthtree);
	
	// print the attributes of the tree before computing "r"
	std::cout << "/----- after building tree " << std::endl; 
	printNode(tree[0]);
	printNode(tree[1]);
	printNode(tree[2]);
	printNode(tree[3]);
	printNode(tree[4]);
	
	// check results 
	std::cout << std::endl << std::endl
			  << "Compare above r to r computed in help-script in R" << std::endl 
			  << "value : r = 0.6348593" << std::endl; 

	return result; 
	
}

int main(){
	
/*	std::cout << "/------------------ test initialization :" << std::endl; 
	testInitialization();

	std::cout << "/------------------ test extent :" << std::endl; 
	testComputeExtent();
	
	std::cout << "/------------------ test morton :" << std::endl; 
	testMorton();
	
	std::cout << "/------------------ test assign particles :" << std::endl; 
	testAssignParticles();
	
	std::cout << "/------------------ test center of mass :" << std::endl; 
	testCenterOfMass();
	
	std::cout << "/------------------ test build :" << std::endl; 
	testBuild();*/
	
	testBuild();
	
	return 0; 
	
}
