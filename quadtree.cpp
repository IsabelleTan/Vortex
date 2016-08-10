//
// Created by Isabelle Tan on 27-07-16.
//

#include <iostream>
#include <cmath>
#include "quadtree.h"

/*
 * help-function for "compute_r"
 */
inline value_type distance_to_COM(const value_type* const xsorted, const value_type* const ysorted, Node const& mynode, const int particle_index)
{
	return sqrt( (xsorted[particle_index] - mynode.xcom)*(xsorted[particle_index] - mynode.xcom) + (ysorted[particle_index] - mynode.ycom)*(ysorted[particle_index] - mynode.ycom) );
} 

/*
 * help-function for "radius": compute r, the radius of the smallest circle centered at the node's COM that contains all particles of that node 
 * given a certain Node
 */
void compute_r(const value_type* const xsorted, const value_type* const ysorted, Node* mynode)
{
	value_type r(distance_to_COM(xsorted, ysorted, *mynode, mynode->part_start)); 
	
	for(size_t i(mynode->part_start+1); i<=mynode->part_end; i++){	
		value_type comp(distance_to_COM(xsorted, ysorted, *mynode, i));
		if(r < comp){r = comp;}
	}
	
	mynode->r = r; 
}

/*
 * help-function for "split"
 * assign radius values to an array of nodes
 * (written like the function centerOfMass() )
 */
void radius(Node* children, value_type* xsorted, value_type* ysorted)
{
	for (int n = 0; n <= 3; ++n) {
		// Set the bounds of the particle indices contained in node n
        int part_start = children[n].part_start;
        int part_end = children[n].part_end;

        // if the node is empty
        if ( (part_start == -1) || (part_end == -1) ){
            // set radius to -1
            children[n].r = -1;			
		}
		// if the node contains only 1 particle
        else if(part_start == part_end) { 
            // set radius to 0
            children[n].r = 0;
		// there is <1 particle in the node
        } else { 
			// compute radius and assign it to child node 				
			value_type r(distance_to_COM(xsorted, ysorted, children[n], children[n].part_start)); 
			for(size_t i(part_start+1); i<=part_end; i++){	
					value_type comp(distance_to_COM(xsorted, ysorted, children[n], i));
					if(r < comp){r = comp;}
			}
			children[n].r = r; 			
		}
	}
}

/*
 * The function that builds the quad-tree
 */
void build(const value_type* const x, const value_type* const y, const value_type* const mass, const int N, const int k, value_type* xsorted, value_type* ysorted, value_type* mass_sorted, Node* tree, int depth){

    // Allocate the index array and compute and store the morton indices in it.
    unsigned int *index = new unsigned int[N];
    double xmin, ymin, ext;
    extent(N, x, y, xmin, ymin, ext);

    morton(N, x, y, xmin, ymin, ext, index, depth);

    // Sort the indices and store the corresponding permutation in keys
    unsigned int *keys = new unsigned int[N];

    // Put the values in keys
    for (int j = 0; j < N; ++j) {
        keys[j] = j;
    }

    sort(N, index, keys);

    // Sort the remaining arrays with the same permutation
    reorder(N, keys, x, y, mass, xsorted, ysorted, mass_sorted);



    // Build the tree

    // Compute the center of mass and the total mass
    value_type xCom = 0, yCom = 0, total_mass = 0;
    for (int i = 0; i < N; ++i) {
        xCom += mass_sorted[i] * xsorted[i];
        yCom += mass_sorted[i] * ysorted[i];
        total_mass += mass_sorted[i];
    }
    xCom /= total_mass;
    yCom /= total_mass;
    
    // leave r initialized to zero.
    // this attribute will never be needed on the first level 
    
    // Allocate variable to contain the index of the next free spot in the nodes array
    int newNodeIndex = 1; // since the root node lies at index 0

    // Create the root node (initially a leaf node)
    tree[0] = Node {0, // root is at level 0
                    0, // morton index
                    -1, // child_id
                    0, // part_start
                    N-1, // part_end
                    total_mass, // node mass
                    xCom, // x center of mass
                    yCom  // y center of mass
    };

    // Subdivide as long as there are more than k particles in the cells
    if (N > k && depth > 0) {
        // There are more than k particles in the node and we haven't reached the maximum depth so split it in four
        // TODO use omp tasking to parallelize this.
        split(tree, tree, depth, index, xsorted, ysorted, mass_sorted, k, &newNodeIndex);
    }
    else {
        // There are less than or equal to k particles in the node
        // or we reached the maximum depth => the node is a leaf
        std::cout << "The root node was not split." << std::endl;
        std::cout << "Particles in root node: " << N << ". Stopping criterion particle number: " << k << std::endl;
        std::cout << "depth = " << depth << std::endl;
    }
    
    delete[] keys;
    delete[] index;
}

/*
 * A function that splits the parent node and creates 4 new children nodes in the array of nodes named "tree".
 */
void split(Node* parent, Node* tree, int depth, unsigned int* index, value_type* xsorted, value_type* ysorted, value_type* mass_sorted, int k, int* newNodeIndex){
// Capture and update the newNodeIndex atomically to avoid race condition
#pragma omp atomic capture
    {
        parent->child_id = *newNodeIndex;
        *newNodeIndex += 4;
    }

    // Compute the level of the children
    int children_level = parent->level +1;

    // Compute the indexvalue of this level so we can easily compute the morton-id's of the children
    int indexValue_level = pow(2,2*(depth - children_level));

    // Initialize the children nodes
    Node child_0 = Node {children_level, // level
                         parent->morton_id, // morton index
                         -1, // child_id
                         -1, // part_start
                         -1, // part_end
                         1, // node mass
                         1, // x center of mass
                         1  // y center of mass
    };

    Node child_1 = Node {children_level, // level
                         parent->morton_id + indexValue_level, // morton index
                         -1, // child_id
                         -1, // part_start
                         -1, // part_end
                         1, // node mass
                         1, // x center of mass
                         1  // y center of mass
    };

    Node child_2 = Node {children_level, // level
                         parent->morton_id + 2*indexValue_level, // morton index
                         -1, // child_id
                         -1, // part_start
                         -1, // part_end
                         1, // node mass
                         1, // x center of mass
                         1  // y center of mass
    };

    Node child_3 = Node {children_level, // level
                         parent->morton_id + 3*indexValue_level, // morton index
                         -1, // child_id
                         -1, // part_start
                         -1, // part_end
                         1, // node mass
                         1, // x center of mass
                         1  // y center of mass
    };

    if(parent->level == 0){
        // Print the maximum index
        unsigned int max = index[0];
        int min = 0;
        for (int l = 0; l < parent->part_end - parent->part_start; ++l) {
            if(index[l]>max){
                max = index[l];
            }
            if(index[l] < min){
                min = index[l];
            }
        }
//        std::cout << "Biggest Morton Index = " << max << std::endl;
//        std::cout << "Smallest Morton Index = " << min << std::endl;
//        std::cout << "Morton limit tree = " <<  child_3.morton_id - 1 + indexValue_level << std::endl;
//        std::cout << "IndexValue at level = " << children_level << " is " << indexValue_level << std::endl;
    }



    // Set a pointer to the first child node
    Node* children = tree + parent->child_id;

    // Put the children nodes into the array at indices [ children[0], ... , children[3] ].
    children[0] = child_0;
    children[1] = child_1;
    children[2] = child_2;
    children[3] = child_3;


    // Assign the particles to the children
    assignParticles(parent, children, depth, index);

    // Assign the total and center of mass to the children, as well as their radius r 
    centerOfMass(children, xsorted, ysorted, mass_sorted);
    radius(children, xsorted, ysorted);

    // TODO Add calls to the p2e() function for each child so that rexps and iexps are computed and set

    // Check number of particles in the children nodes
    for (int i = 0; i < 4; ++i) {
        if (children[i].part_end - children[i].part_start + 1 > k && children[i].level < depth) {
            // There are more than k particles in the node and we haven't reached the maximum depth so split it in four
            split(children + i, tree, depth, index, xsorted, ysorted, mass_sorted, k, newNodeIndex);
        }
        else {
            // There are less than or equal to k particles in the node
            // or we reached the maximum depth => the node is a leaf
        }
    }

}

/*
 * A function that loops over the particles in the parent node and extracts the part_start and part_end values for its
 * child nodes.
 */
void assignParticles(Node* parent, Node* children, int depth, unsigned int* index){
    // TODO the assigment of particles is not bulletproof yet... c sometimes goes beyond 3
    // Check if the parent node is empty
    if(parent->part_start < 0 && parent->part_end <0){
        std::cout << "Something went wrong: the start and end indices of the parents particles are < 0, so we cannot assign any particles to its child nodes." << std::endl;
    }

    // Prepare the variables
    int part_start = parent->part_start;
    int part_end;

    int indexValue_level = pow(2,2*(depth - children[0].level));

    // Start with testing the particles against child node 0
    int c = 0;
    bool cWentOverLimit = false;

    // Loop over all particles
    for (int i = parent->part_start; i <= parent->part_end; ++i) {
//        std::cout << "Level = " << parent->level << ". Checking particle " << i << " with MortonID " << index[i] << " against child node " << c << " with MortonID " << children[c].morton_id << std::endl;

        // Check if the morton index of the particle is smaller than the morton index of the first child node. If so
        // something went wrong.
        if (index[i]<children[c].morton_id){
            std::cout << "Something went wrong." << std::endl;
            std::cout << "Trying to appoint particle " << i << " to Child Node " << c << " at level " << children[0].level <<"." << std::endl;
            std::cout << "Morton domain of child node: " << c << " [" << children[c].morton_id << "," << children[c].morton_id + indexValue_level << "]. Morton index of particle:  " << index[i] <<"." << std::endl;
            break;
        }

        // While the particle is not in the current child node "child", check if the current node is then empty or not.
        // Then, repeat for next child node.
        while(index[i] > children[c].morton_id - 1 + indexValue_level){
            // Check if c is still smaller than 4
            if(c >=4){
                std::cout << "c = " << c << " when trying to assign particle " << i << " with morton id " << index[i] << std::endl;
                std::cout << "This is wrong because c should be in [0,1,2,3] to loop over the child nodes" << std::endl;
                std::cout << "Terminating the assignment of particles." << std::endl;
                cWentOverLimit = true;
                break;
            }

            part_end = i-1;
            // Check whether the current child node is empty
            if(part_end<part_start){
                // This child node is empty so assign negative indices
                children[c].part_start=-1;
                children[c].part_end=-1;

                // Test the particle against the next child node
                c = c+1;
//                std::cout << "Continuing to next children node, node " << c << std::endl;
            } else {
                // The node is not empty so the previous particle was the last particle in the current child node so
                // assign part_start and part_end to the current child node.
                children[c].part_start = part_start;
                children[c].part_end = part_end;

                // Check the same particle again for the next child node
                c=c+1;
//                std::cout << "Continuing to next children node, node " << c << std::endl;

                // Set the start of a new particle group to the current particle
                part_start = i;
            }
        }
        // If c>= 4 then also exit this loop because something is wrong.
        if(cWentOverLimit){
            break;
        }
    }
    // Finished looping over the particles in the parent node
//    std::cout << "Finished looping over particles " << parent->part_start << " to " << parent->part_end << " in parent node" << std::endl;

    // Assign the last particle group to the current child node
    part_end = parent->part_end;
    children[c].part_start = part_start;
    children[c].part_end = part_end;

    // Make the left over child nodes empty nodes
    for (int j = c + 1; j <= 3; j++) {
        children[j].part_start = -1;
        children[j].part_end = -1;
    }

}

void centerOfMass(Node* children, value_type* xsorted, value_type* ysorted, value_type* mass_sorted){
    for (int n = 0; n <= 3; ++n) {
        // Set the bounds of the particle indices contained in node n
        int part_start = children[n].part_start;
        int part_end = children[n].part_end;

        // Check if the node is empty
        if (part_start == -1 || part_end == -1) {
            // Empty node
            children[n].mass = 0;
            children[n].xcom = 0;
            children[n].ycom = 0;
        } else {

            value_type xCenter = 0, yCenter = 0, totalMass = 0;

            // Loop over the particles in node n
            for (int p = part_start; p <= part_end; ++p) {
                xCenter += xsorted[p] * mass_sorted[p];
                yCenter += ysorted[p] * mass_sorted[p];
                totalMass += mass_sorted[p];
            }
            // Divide by the total mass to get the weighted average of the x and y coordinates
            xCenter /= totalMass;
            yCenter /= totalMass;

            // Assign the computed values to the child node
            children[n].xcom = xCenter;
            children[n].ycom = yCenter;
            children[n].mass = totalMass;
        }
    }
}

/*
 * A function that prints the attributes of an object of type Node
 */
void printNode(Node node){
    std::cout<< "\n\nNode at address: " << &node			<< std::endl;
    std::cout<< "level: " 				<< node.level 		<< std::endl;
    std::cout<< "morton_id: " 			<< node.morton_id 	<< std::endl;
    std::cout<< "child_id: " 			<< node.child_id 	<< std::endl;
    std::cout<< "part_start: " 			<< node.part_start 	<< std::endl;
    std::cout<< "part_end: "			<< node.part_end 	<< std::endl;
    std::cout<< "radius r: "			<< node.r		 	<< std::endl;
    std::cout<< "mass: " 				<< node.mass 		<< std::endl;
    std::cout<< "xcom: " 				<< node.xcom 		<< std::endl;
    std::cout<< "ycom: " 				<< node.ycom 		<< std::endl;
}

