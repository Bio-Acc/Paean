#pragma once

#include <cstdint>
#include <numeric>

#include "bin.h"
#include "cuda_util.h"

template <typename T>
struct is_greater_than_one {
    __host__ __device__ bool operator()(const T &x) const {
        return x > 1;
    }
};

template <typename T>
struct min_element {
    __host__ __device__
        T
        operator()(const T &x, const T &y) const {
        return x > y ? y : x;
    }
};

// we discard thrust::gather, just use a simple version instead
template <typename T>
__global__ void gather(uint32_t *indices, T *src, T *dest,
                       uint32_t numOfEntry) {
    uint32_t threadId = blockDim.x * blockIdx.x + threadIdx.x;
    if (threadId < numOfEntry) {
        uint32_t srcId = indices[threadId];
        dest[threadId] = src[srcId];
    }
}

// customized scatter function
template <typename T>
__global__ void scatter_if(T *src, T *dest, uint32_t *indices,
                           uint32_t *flags, uint32_t numOfEntry) {
    uint32_t threadId = blockDim.x * blockIdx.x + threadIdx.x;
    if (threadId < numOfEntry) {
        uint32_t destId = indices[threadId] - 1;
        if (threadId < numOfEntry - 1) {
            if (flags[threadId + 1]) {
                dest[destId] = src[threadId];
            }
        } else {
            dest[destId] = src[threadId];
        }
    }
}

/* these functions is used to sort reads with
 * junctions or without junctions by using
 * cub library.
 */
void cubRadixSortKey(uint64_t *, uint64_t *, uint32_t);
void cubRadixSortInterval(d_Gaps &, d_Gaps &, uint32_t);

// junctions
void cubRadixSortJunction(d_Junctions &, d_Junctions &, uint32_t);
d_Junctions thrustSegmentedScanJunction(d_Junctions &, uint32_t &);

// sort bins and ases by using cub library
void cubRadixSortBin(d_Bins &, d_Bins &, h_Bins &, uint32_t);
void cubRadixSortASE(d_ASEs &, d_ASEs &, uint32_t);

// cub reduce sum
void cubReduceSum(float *, float *, uint32_t);