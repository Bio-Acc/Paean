#include <algorithm>
#include <iostream>
#include <fstream>

#include "fusion.h"

// define bin_read_map
UMAP_PAIR_K bin_read_map;

void Detect_fusion(h_Reads &h_reads, h_Bins &h_bins, 
                   uint32_t numOfRead, uint32_t numOfBin,
                   int max_gap)
{
    for (uint32_t readId = 0; readId < numOfRead; readId++) {
        uint64_t read_s = h_reads.start_[readId];
        uint64_t read_e = h_reads.end_[readId];
        auto its = std::upper_bound(
            h_bins.start_.begin(), h_bins.start_.end(), read_s);
        int idxs = std::distance(h_bins.start_.begin(), its) - 1;
        auto ite = std::upper_bound(
            h_bins.start_.begin(), h_bins.start_.end(), read_e);
        int idxe = std::distance(h_bins.start_.begin(), ite) - 1;
        // three conditions:
        // 1. idxs is not equal to idxe
        // 2. idxs and idxe are both in the bound
        // 3. gap between bins_e and bine_s is greater than max_gap
        if ((idxs != idxe) &&
            (idxs >= 0 && idxs < (int)numOfBin) && 
            (idxe >= 0 && idxe < (int)numOfBin)) {
            if (h_bins.start_[idxe] - h_bins.end_[idxs] > (uint64_t)max_gap) {
                auto key = std::make_pair(h_bins.core[idxs].gid_h, 
                                          h_bins.core[idxe].gid_h);
                bin_read_map[key]++;
            }
        }
    }
}
