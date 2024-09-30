cmake_minimum_required(VERSION 2.8)

find_package(CUDA REQUIRED)

IF(USE_CPP11)
    list(APPEND CUDA_NVCC_FLAGS --gpu-architecture sm_35 -Xptxas -v -use_fast_math -Xptxas -abi=no -std=c++11)
ELSE(USE_CPP11)
    list(APPEND CUDA_NVCC_FLAGS --gpu-architecture sm_35 -Xptxas -v -use_fast_math -Xptxas -abi=no)
ENDIF(USE_CPP11)

set(GPU_CUDA_SRC 
    cudaConstant.cuh
    cudaReductions.cuh
    customLDG.cuh 
    fakeComplex.cuh 
    gFunctions.cuh 
    globals.h
    gpu.cuh 
    HC_SF_GPU.cu 
#   HC_SF_GPU.h 
    memoryAlloc2.cuh 
    realSphericalHarmonics.cuh
    types.h 
)

CUDA_ADD_LIBRARY(discamb_cuda ${GPU_CUDA_SRC})


