#ifndef PAEAN_SORT_READ_CUH
#define PAEAN_SORT_READ_CUH

#include <cstdint>
#include <numeric>

#include <thrust/sort.h>
#include <thrust/scan.h>
#include <thrust/fill.h>
#include <thrust/copy.h>
#include <thrust/reduce.h>
#include <thrust/replace.h>
#include <thrust/gather.h>
#include <thrust/scatter.h>
#include <thrust/host_vector.h>

#include "bin.h"

template <typename T>
struct is_greater_than_one {
    bool operator()(const T &x) const {
        return x > 1;
    }
};

template <typename T>
struct min_element {
    T operator()(const T &x, const T &y) const {
        return x > y ? y : x;
    }
};

// customized scatter function
template <typename T>
void scatter_if(T *src, T *dest, uint32_t *indices,
                uint32_t *flags, uint32_t numOfEntry) {
#pragma omp parallel for
    for (uint32_t i = 0; i < numOfEntry; i++) {
        uint32_t destId = indices[i] - 1;
        if (i < numOfEntry - 1) {
            if (flags[i + 1]) {
                dest[destId] = src[i];
            }
        } else {
            dest[destId] = src[i];
        }
    }
}

// junctions
void thrustSortJunction(h_Junctions &);
h_Junctions thrustSegmentedScanJunction(h_Junctions &);

// sort bins and ases by using thrust library
void thrustSortBin(h_Bins &);
void thrustSortASE(h_ASEs &);

#endif // PAEAN_SORT_READ_CUH