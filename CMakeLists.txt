cmake_minimum_required(VERSION 3.5)
project(Vortex)

set(CXX g++-5)                                                                  # Use g++ version 5 (required to use openMP)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp -mavx -mfma -O3")   # Compile with these flags


set(SOURCE_FILES
        main.cpp
        expansion.cpp
        solver.cpp
        quadtree.cpp
        solver_test.cpp
        fields.cpp
        fields_test.cpp
        ADIdiffusion.cpp
        ADIdiffusion_test.cpp
        )

set(HEADER_FILES
        expansion.h
        solver.h
        quadtree.h
        solver_test.h
        fields.h
        fields_test.h
        datapoints.h
        ADIdiffusion.h
        ADIdiffusion_test.h
        )

add_executable(VORTEX ${SOURCE_FILES} ${HEADER_FILES})

##################
### Compile flags
##################

# -std=c++11                                # Use c++11

# -fopenmp                                  # Use OpenMP

# -ftree-vectorize -fopt-info-vec -O        # Tell the compiler to vectorize code if it thinks it's possible,
                                            # and show output about loops that are vectorized. The -O is necessary.

# -mavx -mfma                               # Use AVX intrinsics and use Fused Multiply Add (FMA) intrinsics for
                                            # vectorized code.