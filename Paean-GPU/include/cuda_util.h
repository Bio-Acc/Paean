#pragma once

#include <cstdio>
#include <cstdlib>

#include "cuda_runtime.h"

#define DIV_UP(x, y) ((x + y - 1) / y)

#define CUDA_SAFE_CALL(call)                                              \
    {                                                                     \
        cudaError_t err = call;                                           \
        if (cudaSuccess != err) {                                         \
            fprintf(stderr, "Cuda error in file '%s' in line %i : %s.\n", \
                    __FILE__, __LINE__, cudaGetErrorString(err));         \
            exit(EXIT_FAILURE);                                           \
        }                                                                 \
    }
