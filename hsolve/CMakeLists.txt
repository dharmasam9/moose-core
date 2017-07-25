cmake_minimum_required(VERSION 2.6)
include_directories(../basecode ../utility ../kinetics ../external/debug)

if(WITH_CUDA)
    add_definitions(-DUSE_CUDA)
    #SET(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS};-gencode arch=compute_30,code=sm_30)
    #SET(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS};-std=c++11 )
    #SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC" )
    SET(CMAKE_CXX_FLAGS "-fPIC" )
    FILE(GLOB CUDA_SRC *.cu)
    INCLUDE_DIRECTORIES(${CUDA_INCLUDE_DIRS})
    CUDA_COMPILE(HSOLVE_CUDA_OBJS ${CUDA_SRC})
    CUDA_ADD_LIBRARY(hsolve_cuda ${HSOLVE_CUDA_OBJS})
    TARGET_LINK_LIBRARIES(hsolve_cuda ${CUDA_LIBRARIES} ${CUDA_cusparse_LIBRARY} ${CUDA_cublas_LIBRARY} )
endif()

add_library(hsolve
    	HSolveStruct.cpp
	HinesMatrix.cpp
	HSolvePassive.cpp
	RateLookup.cpp
	HSolveActive.cpp
	HSolveActiveSetup.cpp
	HSolveInterface.cpp
	HSolve.cpp
	HSolveUtils.cpp
	testHSolve.cpp
	ZombieCompartment.cpp
	ZombieCaConc.cpp
	ZombieHHChannel.cpp
)

if(WITH_CUDA)
    TARGET_LINK_LIBRARIES(hsolve hsolve_cuda)
endif(WITH_CUDA)

