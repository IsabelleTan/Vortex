//
// Created by Isabelle Tan on 04-05-16.
//

#include <iostream>
#include "morton.h"
#include "quadtree.h"
#include "expansion.h"

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

// Test the morton indices against a test set of 7 particles with known indices.
bool testMorton(){
    int N = 7;
    int depth = 2;
    // Set the resulting boolean to true
    bool result = true;

    // Initiate arrays
    // Create particle coordinate arrays
    value_type x [N] = {0.1, 0.4, 0.8, 0.9, 0.9, 0, 1};
    value_type y [N] = {0.3, 0.6, 0.4, 0.8, 0.1, 0, 1};
    int *index = new int[N];
    uint32_t control [N] = {8, 3, 13, 5, 15, 10, 5};
    double xmin, ymin, ext;

    // Compute indices
    morton(N, x, y, xmin, ymin, ext, index);

    // test the resulting indices against the control indices
    for (int i = 0; i < N; ++i) {

        if(index[i] != control[i]){
            result = false;
            std::cout << "Test failed for particle " << i << std::endl;
        }
    }

    // Print the results
    std::cout<< "mortonSerial: \t" << " control"<< std::endl;
    for (int j = 0; j < N; ++j) {
        std::cout << index[j] << "\t\t\t\t\t" << control[j] << std::endl;
    }

    // Free the memory
    delete[] index;

    // Print in case of success.
    if (result) {
        std::cout << "Test succeeded." << std::endl;
    }

    return result;
}

// Test the assignment of part_start and part_end to child nodes
bool testAssignParticles(){
    // Don't change these variables.
    int N = 7;
    int depth = 2;
    // Set the resulting boolean to true
    bool result = true;

    // Initiate arrays
    // Create index array
    int index_0[N] = {3, 5, 5, 8, 10, 13, 15};
    int index_1[N] = {4, 5, 5, 8, 9, 10, 11};
    int index_2[N] = {3, 13, 13, 13, 14, 15, 15};

    // Make arrays containing the known true values
    int part_start_control_0[4] = {0, 1, 3, 5};
    int part_end_control_0[4] = {0, 2, 4, 6};
    int part_start_control_1[4] = {-1, 0, 3, -1};
    int part_end_control_1[4] = {-1, 2, 6, -1};
    int part_start_control_2[4] = {0, -1, -1, 1};
    int part_end_control_2[4] = {0, -1, -1, 6};

    // Make a root node
    Node root {0, // level
               0, // morton index
               0, // child_id
               0, // part_start
               6, // part_end
               0, // node mass
               0, // x center of mass
               0  // y center of mass
    };


    Node children[4];

    // Compute the indexvalue of this level so we can easily compute the morton-id's of the children
    int indexValue_level = pow(2,2*(depth - (root.level + 1)));

    // Initialize the children nodes
    Node child_0 = Node {root.level + 1, // level
                         root.morton_id, // morton index
                         -1, // child_id
                         -1, // part_start
                         -1, // part_end
                         1, // node mass
                         1, // x center of mass
                         1  // y center of mass
    };

    Node child_1 = Node {root.level + 1, // level
                         root.morton_id + indexValue_level, // morton index
                         -1, // child_id
                         -1, // part_start
                         -1, // part_end
                         1, // node mass
                         1, // x center of mass
                         1  // y center of mass
    };

    Node child_2 = Node {root.level + 1, // level
                         root.morton_id + 2 * indexValue_level, // morton index
                         -1, // child_id
                         -1, // part_start
                         -1, // part_end
                         1, // node mass
                         1, // x center of mass
                         1  // y center of mass
    };

    Node child_3 = Node {root.level + 1, // level
                         root.morton_id + 3 * indexValue_level, // morton index
                         -1, // child_id
                         -1, // part_start
                         -1, // part_end
                         1, // node mass
                         1, // x center of mass
                         1  // y center of mass
    };

    // Put the children nodes into the array
    children[0] = child_0;
    children[1] = child_1;
    children[2] = child_2;
    children[3] = child_3;


    // Test set 0
    assignParticles(&root, children, depth, index_0);

    // Check the start and end indices of this node
    for (int i = 0; i < 4; ++i) {
        if ((children[i].part_start != part_start_control_0[i]) || (children[i].part_end != part_end_control_0[i])) {
            result = false;
            std::cout << "Test failed for set 0 and node " << i << std::endl;
        }
    }

    // Test set 1
    assignParticles(&root, children, depth, index_1);

    // Check the start and end indices of this node
    for (int i = 0; i < 4; ++i) {
        if ((children[i].part_start != part_start_control_1[i]) || (children[i].part_end != part_end_control_1[i])) {
            result = false;
            std::cout << "Test failed for set 1 and node " << i << std::endl;
            std::cout << "children[i].part_start = " << children[i].part_start << ". children[i].part_end = " << children[i].part_end << std::endl;
            std::cout << "node " << i << " morton_id = " << children[i].morton_id << " non-inclusive upper bound = " << children[i].morton_id + indexValue_level << std::endl;
        }
    }

    // Test set 2
    assignParticles(&root, children, depth, index_2);

    // Check the start and end indices of this node
    for (int i = 0; i < 4; ++i) {
        if ((children[i].part_start != part_start_control_2[i]) || (children[i].part_end != part_end_control_2[i])) {
            result = false;
            std::cout << "Test failed for set 2 and node " << i << std::endl;
        }
    }

    // Print in case of success.
    if (result) {
        std::cout << "Test succeeded." << std::endl;
    }

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
    int depth = 2;
    int k = 8;

    value_type x[7] = {0.1, 0.15, 0.4, 0.6, 0.7, 0.71, 0.8};
    value_type y[7] = {0.6, 0.85, 0.65, 0.3, 0.6, 0.1, 0.4};
    value_type mass[7] = {1, 2, 3, 4, 5, 6, 7};

    value_type xsorted[7];
    value_type ysorted[7];
    value_type mass_sorted[7];

    // Allocate the tree array containing all the nodes
    int maxNodes = (int) std::min((float)8 * N / k, (float)pow(4, depth));
    Node* tree = new Node[maxNodes];

    build(x, y, mass, N, k, xsorted, ysorted, mass_sorted, tree, depth);

    // Test the resulting tree
    if (tree[0].child_id != -1) {
        std::cout << "Test of stopping criterion failed, the root node is not a leaf node." << std::endl;
        std::cout << "root.child_id = " << tree[0].child_id << std::endl;
        result = false;
    } else {
        std::cout << "Root node is an empty node => The test for the stopping criterion passed." << std::endl;
    }


    // Test building the tree for k=2
    std::cout << "Testing the build() function for a tree with N=7, k=2 and depth=2." << std::endl;
    k = 2;
    build(x, y, mass, N, k, xsorted, ysorted, mass_sorted, tree, depth);

    // Test the resulting tree
    for (int i = 0; i < 13; ++i) {
        printNode(tree[i]);
    }

    if (result) {
        std::cout << "Check the created nodes to validate the function." << std::endl;
    }

    delete[] tree;
    return result;
}

// Test the computation of attribute "r"
bool testr(){
	
	bool result = true; 
	
	// Create particles and particle arrays 
	int N = 10;
    int depth = 2;
    int k = 8;

    value_type x[10] = {0.1, 0.15, 0.4, 0.6, 0.7, 0.71, 0.8, 0.74, 0.48, 0.41};
    value_type y[10] = {0.6, 0.85, 0.65, 0.3, 0.6, 0.1, 0.4, 0.33, 0.42, 0.57};
    value_type mass[10] = {1, 2, 3, 4, 5, 6, 7, 2.3, 1.4, 5.7};

    value_type xsorted[10];
    value_type ysorted[10];
    value_type mass_sorted[10];

    // Allocate the tree array containing all the nodes
    int maxNodes = (int) std::min((float)8 * N / k, (float)pow(4, depth));
    Node* tree = new Node[maxNodes];

	// create a tree to be tested 
    build(x, y, mass, N, k, xsorted, ysorted, mass_sorted, tree, depth);
	
	// print the attributes of the tree before computing "r"
	std::cout << "/----- after building tree " << std::endl; 
	printNode(tree[0]);
	printNode(tree[1]);
	printNode(tree[2]);
	printNode(tree[3]);
	printNode(tree[4]);
	printNode(tree[5]);
	printNode(tree[6]);
	printNode(tree[7]);
	printNode(tree[8]);
	printNode(tree[9]);
	
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
	
	testr();
	
	return 0; 
	
}
