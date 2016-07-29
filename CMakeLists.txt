cmake_minimum_required(VERSION 3.5)
project(Vortex)

set(CXX g++-5)                                                          # Use g++ version 5 (required to use openMP)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp -pthread -mavx -mfma -O3")   # Compile with these flags


set(SOURCE_FILES
        main.cpp
        simulation.cpp
        datapoints.cpp
        datapoints_test.cpp
        expansion.cpp
        expansion_test.cpp
        solver.cpp
        morton.cpp
        morton_test.cpp
        morton_time.cpp
        quadtree.cpp
        quadtree_test.cpp
        quadtree_time.cpp
        multipole.cpp
        multipole_test.cpp
        solver_test.cpp
        fields.cpp
        fields_test.cpp
        )

set(HEADER_FILES
		simulation.h
		datapoints.h
        expansion.h
        multipole.h
        solver.h
        morton.h
        quadtree.h
        fields.h
        )

add_executable(datapoints_test datapoints_test.cpp datapoints.cpp)		# tested: OK 
add_executable(morton_time morton_time.cpp morton.cpp datapoints.cpp) 	# tested: OK 
add_executable(morton_test morton_test.cpp morton.cpp datapoints.cpp) 	# tested: OK
add_executable(quadtree_test quadtree_test.cpp morton.cpp quadtree.cpp) # tested: OK, except for very weird error mesage all the way at the end of execution ... ???
add_executable(quadtree_time quadtree_time.cpp morton.cpp quadtree.cpp) # tested: OK
add_executable(expansion_test expansion_test.cpp expansion.cpp morton.cpp complexAVX.cpp) 		# tested: OK
add_executable(multipole_test multipole_test.cpp multipole.cpp expansion.cpp morton.cpp quadtree.cpp complexAVX.cpp) 

#add_executable(VORTEX ${SOURCE_FILES} ${HEADER_FILES}) 				# Shoshi: I'm sorry, I commented this because it's a bit mysterious to me ... I prefer writing each build command separately so I can compile one executable at a time if needed, and I have a better overview of the dependancies ... 
 
#---------------------#
#   Compile flags     #
#---------------------#

# -std=c++11                                # Use c++11

# -pthread									# use multithreading (for extent() in quadtree)

# -fopenmp                                  # Use OpenMP

# -ftree-vectorize -fopt-info-vec -O        # Tell the compiler to vectorize code if it thinks it's possible,
                                            # and show output about loops that are vectorized. The -O is necessary.

# -mavx -mfma                               # Use AVX intrinsics and use Fused Multiply Add (FMA) intrinsics for
                                            # vectorized code.
