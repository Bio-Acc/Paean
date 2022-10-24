#ifndef PAEAN_BIN_ASSIGN_KERNEL_H
#define PAEAN_BIN_ASSIGN_KERNEL_H

#include "bin.h"
#include "bisect.h"

/* this kernel is used to find and count the reads
 * which overlap genes.
 */
void assign_read_kernel(h_Bins &h_bins, h_Reads &h_reads,
                        h_Gaps &h_gap_singles, h_Gaps &h_gap_pairs,
                        uint32_t *h_readCount, bool paired_end) {
    uint32_t numOfBin = h_bins.start_.size();
    uint32_t numOfRead = h_reads.start_.size();
    uint32_t numOfGap = h_gap_singles.start_.size();

#pragma omp parallel for
    for (uint32_t binId = 0; binId < numOfBin; binId++) {
        int overlap = count_overlap(h_bins.start_[binId],
                                    h_bins.end_[binId],
                                    h_reads.start_.data(),
                                    h_reads.end_.data(),
                                    numOfRead);
        uint32_t gap = 0;
        if (paired_end) {
            gap = count_outer_contain(h_bins.start_[binId], h_bins.end_[binId],
                                      h_gap_singles.start_.data(),
                                      h_gap_singles.end_.data(),
                                      numOfGap);
            // remove duplicates
            gap += count_inner_contain(h_bins.start_[binId], h_bins.end_[binId],
                                       h_gap_pairs.start_.data(),
                                       h_gap_pairs.end_.data(),
                                       numOfGap);
        }
        uint32_t readCount = overlap - gap;
        h_readCount[binId] = readCount;
        h_bins.core[binId].readCount += readCount;
        // #define DEBUG
#ifdef DEBUG
        printf("read count: %u\n", h_bins.core[binId].readCount);
#endif
    }
}

/* this kernel is used to find and count the ases
 * for each gene.
 */
void assign_ASE_kernel(h_Bins &h_bins, h_ASEs &h_ases, h_Assist &h_assist) {
    uint32_t numOfBin = h_bins.start_.size();
    uint32_t numOfASE = h_ases.start_.size();

#pragma omp parallel for
    for (uint32_t binId = 0; binId < numOfBin; binId++) {
        find_inner_contain(h_bins.start_[binId], h_bins.end_[binId],
                           h_ases.start_.data(), numOfASE,
                           h_assist.start_.data(), h_assist.end_.data(), binId);
        for (uint32_t aseId = h_assist.start_[binId]; aseId < h_assist.end_[binId]; aseId++) {
            if (h_ases.end_[aseId] <= h_bins.end_[binId]) {
#pragma omp critical
                {
                    uint32_t binCount = h_ases.core[aseId].bin_h.binCount;
                    // only store up to `binNameSize` gene names
                    if (binCount < binNameSize) {
                        h_ases.core[aseId].bin_h.bins[binCount] = h_bins.core[binId].gid_h;
                        h_ases.core[aseId].bin_h.binCount++;
                    }
                }
            }
        }
    }
}

/* we use junction table and binary search algorithm
 * to assign reads' junctions to ases.
 */
void assign_read_ASE_kernel(h_ASEs h_ases, h_Junctions h_junctions,
                            ASEPsi *h_ase_psi) {
    uint32_t numOfASE = h_ases.start_.size();
    uint32_t numOfJunction = h_junctions.start_.size();

#pragma omp parallel for
    for (uint32_t aseId = 0; aseId < numOfASE; aseId++) {
        float countIn = 0, countOut = 0;
        int idx;

        // out junction
        for (uint32_t i = 0; i < h_ases.core[aseId].coordinateCountOut; i += 2) {

            idx = find_coincide(h_ases.core[aseId].coordinates[i],
                                h_junctions.start_.data(), numOfJunction);
            if (idx != -1) {
                if (h_ases.core[aseId].coordinates[i + 1] == h_junctions.end_[idx]) {
                    countOut += h_junctions.count[idx];
                } else {
                    // continue to search
                    while (++idx < numOfJunction) {
                        if (h_ases.core[aseId].coordinates[i] != h_junctions.start_[idx]) {
                            break;
                        }
                        if (h_ases.core[aseId].coordinates[i + 1] == h_junctions.end_[idx]) {
                            countOut += h_junctions.count[idx];
                            break;
                        }
                    }
                }
            }
        }

        // in junctions, skip out junction
        for (uint32_t j = 0; j < h_ases.core[aseId].coordinateCountIn; j += 2) {

            uint32_t k = j + h_ases.core[aseId].coordinateCountOut;
            idx = find_coincide(h_ases.core[aseId].coordinates[k],
                                h_junctions.start_.data(), numOfJunction);
            if (idx != -1) {
                if (h_ases.core[aseId].coordinates[k + 1] == h_junctions.end_[idx]) {
                    countIn += h_junctions.count[idx];
                } else {
                    // continue to search
                    while (++idx < numOfJunction) {
                        if (h_ases.core[aseId].coordinates[k] != h_junctions.start_[idx]) {
                            break;
                        }
                        if (h_ases.core[aseId].coordinates[k + 1] == h_junctions.end_[idx]) {
                            countIn += h_junctions.count[idx];
                            break;
                        }
                    }
                }
            }
        }
        // store into d_ase_psi
        h_ase_psi[aseId] = ASEPsi{h_ases.core[aseId].gid_h,
                                  h_ases.core[aseId].bin_h,
                                  countIn,
                                  countOut,
                                  0,
                                  0,
                                  0};
    }
}

// only designed for RI event
void assign_read_RI_kernel(h_ASEs h_ases, h_Reads h_nj_reads,
                           h_Junctions h_junctions, ASEPsi *h_ase_psi) {
    uint32_t numOfASE = h_ases.start_.size();
    uint32_t numOf_nj_Read = h_nj_reads.start_.size();
    uint32_t numOfJunction = h_junctions.start_.size();

#pragma omp parallel for
    for (uint32_t aseId = 0; aseId < numOfASE; aseId++) {
        float countIn = 0, countOut = 0;
        int idx;

        for (uint32_t i = 0; i < h_ases.core[aseId].coordinateCountOut; i += 2) {

            idx = find_coincide(h_ases.core[aseId].coordinates[i],
                                h_junctions.start_.data(), numOfJunction);
            if (idx != -1) {
                if (h_ases.core[aseId].coordinates[i + 1] == h_junctions.end_[idx]) {
                    countOut += h_junctions.count[idx];
                } else {
                    // continue to search
                    while (++idx < numOfJunction) {
                        if (h_ases.core[aseId].coordinates[i] != h_junctions.start_[idx]) {
                            break;
                        }
                        if (h_ases.core[aseId].coordinates[i + 1] == h_junctions.end_[idx]) {
                            countOut += h_junctions.count[idx];
                            break;
                        }
                    }
                }
            }
        }

        // in junctions, skip out junction
        for (uint32_t j = 0; j < h_ases.core[aseId].coordinateCountIn; j += 2) {

            uint32_t k = j + h_ases.core[aseId].coordinateCountOut;

            int nj_overlap = count_overlap(h_ases.core[aseId].coordinates[k],
                                           h_ases.core[aseId].coordinates[k + 1],
                                           h_nj_reads.start_.data(),
                                           h_nj_reads.end_.data(),
                                           numOf_nj_Read);
            countIn += nj_overlap;
        }

        countIn /= 2;

        h_ase_psi[aseId] = ASEPsi{h_ases.core[aseId].gid_h,
                                  h_ases.core[aseId].bin_h,
                                  countIn,
                                  countOut,
                                  0,
                                  0,
                                  0};
    }
}

#endif // PAEAN_BIN_ASSIGN_KERNEL_H
