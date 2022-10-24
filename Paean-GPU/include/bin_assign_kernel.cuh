#pragma once

#include "bin.h"
#include "bisect.cuh"

#include <stdio.h>
#include <assert.h>

/* this kernel is used to find and count the reads
 * which overlap genes, including reads with junctions
 * and reads without junctions.
 */
__global__ void gpu_assign_read_kernel(d_Bins d_bins, uint32_t numOfBin,
                                       d_Reads d_reads, uint32_t numOfRead,
                                       d_Gaps d_gap_singles, d_Gaps d_gap_pairs,
                                       uint32_t numOfGap, uint32_t *d_readCount,
                                       bool paired_end) {
    uint32_t binId = blockDim.x * blockIdx.x + threadIdx.x;

    if (binId < numOfBin) {
        int overlap = count_overlap(d_bins.start_[binId], d_bins.end_[binId],
                                    d_reads.start_, d_reads.end_,
                                    numOfRead);
        uint32_t gap = 0;
        if (paired_end) {
            gap = count_outer_contain(d_bins.start_[binId], d_bins.end_[binId],
                                      d_gap_singles.start_, d_gap_singles.end_,
                                      numOfGap);
            // remove duplicates
            gap += count_inner_contain(d_bins.start_[binId], d_bins.end_[binId],
                                       d_gap_pairs.start_, d_gap_pairs.end_,
                                       numOfGap);
        }
        // Note: these are partial assignments, we need to merge
        // d_reads and d_nj_reads to the complete d_reads
        uint32_t readCount = overlap - gap;
        d_readCount[binId] = readCount;
        d_bins.core[binId].readCount += readCount;
        // #define DEBUG
#ifdef DEBUG
        printf("read count: %u\n", d_bins.core[binId].readCount);
#endif
    }
}

/* this kernel is used to find and count the ases
 * for each gene.
 */
__global__ void gpu_assign_ASE_kernel(d_Bins d_bins, uint32_t numOfBin,
                                      d_ASEs d_ases, uint32_t numOfASE,
                                      Assist d_assist) {
    uint32_t binId = blockDim.x * blockIdx.x + threadIdx.x;

    if (binId < numOfBin) {
        find_inner_contain(d_bins.start_[binId], d_bins.end_[binId],
                           d_ases.start_, numOfASE,
                           d_assist.start_, d_assist.end_, binId);
        for (uint32_t aseId = d_assist.start_[binId];
             aseId < d_assist.end_[binId]; aseId++) {
            if (d_ases.end_[aseId] <= d_bins.end_[binId]) {
                uint32_t binCount = atomicAdd(&(d_ases.core[aseId].bin_h.binCount), 1);
                // only store up to `binNameSize` gene names
                if (binCount < binNameSize) {
                    d_ases.core[aseId].bin_h.bins[binCount] = d_bins.core[binId].gid_h;
                }
            }
        }
    }
}

/* we use junction table and binary search algorithm
 * to assign reads' junctions to ases.
 */
__global__ void gpu_assign_read_ASE_kernel(d_ASEs d_ases, uint32_t numOfASE,
                                           d_Junctions d_junctions, uint32_t numOfJunction,
                                           ASEPsi *d_ase_psi) {
    uint32_t aseId = blockDim.x * blockIdx.x + threadIdx.x;
    float countIn = 0, countOut = 0;
    int idx;

    if (aseId < numOfASE) {
        // out junction
        for (uint32_t i = 0; i < d_ases.core[aseId].coordinateCountOut; i += 2) {

            idx = find_coincide(d_ases.core[aseId].coordinates[i],
                                d_junctions.start_, numOfJunction);
            if (idx != -1) {
                if (d_ases.core[aseId].coordinates[i + 1] == d_junctions.end_[idx]) {
                    countOut += d_junctions.count[idx];
                } else {
                    // continue to search
                    while (++idx < numOfJunction) {
                        if (d_ases.core[aseId].coordinates[i] != d_junctions.start_[idx]) {
                            break;
                        }
                        if (d_ases.core[aseId].coordinates[i + 1] == d_junctions.end_[idx]) {
                            countOut += d_junctions.count[idx];
                            break;
                        }
                    }
                }
            }
        }

        // in junctions, skip out junction
        for (uint32_t j = 0; j < d_ases.core[aseId].coordinateCountIn; j += 2) {

            uint32_t k = j + d_ases.core[aseId].coordinateCountOut;
            idx = find_coincide(d_ases.core[aseId].coordinates[k],
                                d_junctions.start_, numOfJunction);
            if (idx != -1) {
                if (d_ases.core[aseId].coordinates[k + 1] == d_junctions.end_[idx]) {
                    countIn += d_junctions.count[idx];
                } else {
                    // continue to search
                    while (++idx < numOfJunction) {
                        if (d_ases.core[aseId].coordinates[k] != d_junctions.start_[idx]) {
                            break;
                        }
                        if (d_ases.core[aseId].coordinates[k + 1] == d_junctions.end_[idx]) {
                            countIn += d_junctions.count[idx];
                            break;
                        }
                    }
                }
            }
        }

        // store into d_ase_psi
        d_ase_psi[aseId] = ASEPsi{d_ases.core[aseId].gid_h,
                                  d_ases.core[aseId].bin_h,
                                  countIn,
                                  countOut,
                                  0,
                                  0,
                                  0};
    }
}

// only designed for RI event
__global__ void gpu_assign_read_RI_kernel(d_ASEs d_ases, uint32_t numOfASE,
                                          d_Reads d_nj_reads, uint32_t numOf_nj_Read,
                                          d_Junctions d_junctions, uint32_t numOfJunction,
                                          ASEPsi *d_ase_psi) {
    uint32_t aseId = blockDim.x * blockIdx.x + threadIdx.x;
    float countIn = 0, countOut = 0;
    int idx;

    if (aseId < numOfASE) {
        for (uint32_t i = 0; i < d_ases.core[aseId].coordinateCountOut; i += 2) {

            idx = find_coincide(d_ases.core[aseId].coordinates[i],
                                d_junctions.start_, numOfJunction);
            if (idx != -1) {
                if (d_ases.core[aseId].coordinates[i + 1] == d_junctions.end_[idx]) {
                    countOut += d_junctions.count[idx];
                } else {
                    // continue to search
                    while (++idx < numOfJunction) {
                        if (d_ases.core[aseId].coordinates[i] != d_junctions.start_[idx]) {
                            break;
                        }
                        if (d_ases.core[aseId].coordinates[i + 1] == d_junctions.end_[idx]) {
                            countOut += d_junctions.count[idx];
                            break;
                        }
                    }
                }
            }
        }

        // in junctions, skip out junction
        for (uint32_t j = 0; j < d_ases.core[aseId].coordinateCountIn; j += 2) {

            uint32_t k = j + d_ases.core[aseId].coordinateCountOut;

            int nj_overlap = count_overlap(d_ases.core[aseId].coordinates[k],
                                           d_ases.core[aseId].coordinates[k + 1],
                                           d_nj_reads.start_, d_nj_reads.end_,
                                           numOf_nj_Read);
            countIn += nj_overlap;
        }

        countIn /= 2;

        d_ase_psi[aseId] = ASEPsi{d_ases.core[aseId].gid_h,
                                  d_ases.core[aseId].bin_h,
                                  countIn,
                                  countOut,
                                  0,
                                  0,
                                  0};
    }
}
