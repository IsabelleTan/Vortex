cmake_minimum_required(VERSION 3.5)
project(Vortex)

set(CXX g++-5)                                                          # Use g++ version 5 (required to use openMP)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp -pthread -mavx -mfma -O3")   # Compile with these flags


set(SOURCE_FILES
        main.cpp
        datapoints.cpp
        datapoints_test.cpp
        expansion.cpp
        solver.cpp
        quadtree.cpp
        solver_test.cpp
        fields.cpp
        fields_test.cpp
        )

set(HEADER_FILES
		datapoints.h
        expansion.h
        solver.h
        quadtree.h
        solver_test.h
        fields.h
        fields_test.h
        )

#add_executable(VORTEX ${SOURCE_FILES} ${HEADER_FILES}) 				# Shoshi: I'm sorry, I commented this because it's a bit mysterious to me ... I prefer writing each build command separately so I know exactly what it does 
add_executable(datapoints_test datapoints_test.cpp datapoints.cpp)		# tested: OK 
 
##################
### Compile flags
##################

# -std=c++11                                # Use c++11

# -pthread									# use multithreading (for extent() in quadtree)

# -fopenmp                                  # Use OpenMP

# -ftree-vectorize -fopt-info-vec -O        # Tell the compiler to vectorize code if it thinks it's possible,
                                            # and show output about loops that are vectorized. The -O is necessary.

# -mavx -mfma                               # Use AVX intrinsics and use Fused Multiply Add (FMA) intrinsics for
                                            # vectorized code.
