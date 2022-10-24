#include "util.h"
#include "sort.h"

void thrustSortJunction(h_Junctions &h_junctions) {
    thrust::stable_sort_by_key(h_junctions.end_.begin(), h_junctions.end_.end(),
                               h_junctions.start_.begin());
    thrust::stable_sort_by_key(h_junctions.start_.begin(), h_junctions.start_.end(),
                               h_junctions.end_.begin());
}

h_Junctions thrustSegmentedScanJunction(h_Junctions &h_junctions) {
    // segmented prefix sum
    uint32_t numOfJunction = h_junctions.start_.size();
    thrust::host_vector<uint32_t> counts_start(numOfJunction);
    thrust::host_vector<uint32_t> counts_end(numOfJunction);
    thrust::fill(counts_start.begin(), counts_start.end(), 1);
    thrust::fill(counts_end.begin(), counts_end.end(), 1);

    thrust::inclusive_scan_by_key(h_junctions.start_.begin(),
                                  h_junctions.start_.end(),
                                  counts_start.begin(),
                                  counts_start.begin());
    thrust::inclusive_scan_by_key(h_junctions.end_.begin(),
                                  h_junctions.end_.end(),
                                  counts_end.begin(),
                                  counts_end.begin());

    thrust::host_vector<uint32_t> counts(numOfJunction);
    thrust::transform(counts_start.begin(), counts_start.end(),
                      counts_end.begin(), counts.begin(),
                      min_element<uint32_t>());

    // set flags
    thrust::host_vector<uint32_t> flags(numOfJunction);
    thrust::replace_copy_if(counts.begin(), counts.end(), flags.begin(),
                            is_greater_than_one<uint32_t>(), 0);

    // compute offsets
    thrust::host_vector<uint32_t> indices(numOfJunction);
    thrust::inclusive_scan(flags.begin(), flags.end(), indices.begin());

    // calculate new numOfJunction
    uint32_t new_numOfJunction = indices.back();

    h_Junctions h_junctions_out;
    h_junctions_out.start_.resize(new_numOfJunction);
    h_junctions_out.end_.resize(new_numOfJunction);
    h_junctions_out.count.resize(new_numOfJunction);

    scatter_if<uint64_t>(h_junctions.start_.data(),
                         h_junctions_out.start_.data(),
                         indices.data(), flags.data(),
                         numOfJunction);
    scatter_if<uint64_t>(h_junctions.end_.data(),
                         h_junctions_out.end_.data(),
                         indices.data(), flags.data(),
                         numOfJunction);
    scatter_if<uint32_t>(counts.data(), h_junctions_out.count.data(),
                         indices.data(), flags.data(),
                         numOfJunction);

    return h_junctions_out;
}

void thrustSortBin(h_Bins &h_bins) {
    h_Bins h_bins_out;
    uint32_t numOfBin = h_bins.start_.size();
    h_bins_out.start_.resize(numOfBin);
    h_bins_out.end_.resize(numOfBin);
    h_bins_out.strand.resize(numOfBin);
    h_bins_out.core.resize(numOfBin);

    std::vector<uint32_t> indices(h_bins.start_.size());
    std::iota(indices.begin(), indices.end(), 0);

    // sort
    thrust::stable_sort_by_key(h_bins.start_.begin(), h_bins.start_.end(),
                               indices.begin());

    // gather the remaining fields
    thrust::gather(indices.begin(), indices.end(), h_bins.end_.begin(),
                   h_bins_out.end_.begin());
    thrust::gather(indices.begin(), indices.end(), h_bins.strand.begin(),
                   h_bins_out.strand.begin());
    thrust::gather(indices.begin(), indices.end(), h_bins.core.begin(),
                   h_bins_out.core.begin());
    // swap
    h_bins.end_.swap(h_bins_out.end_);
    h_bins.strand.swap(h_bins_out.strand);
    h_bins.core.swap(h_bins_out.core);
}

void thrustSortASE(h_ASEs &h_ases) {
    h_ASEs h_ases_out;
    uint32_t numOfASE = h_ases.start_.size();
    h_ases_out.start_.resize(numOfASE);
    h_ases_out.end_.resize(numOfASE);
    h_ases_out.strand.resize(numOfASE);
    h_ases_out.core.resize(numOfASE);

    std::vector<uint32_t> indices(h_ases.start_.size());
    std::iota(indices.begin(), indices.end(), 0);

    // sort
    thrust::stable_sort_by_key(h_ases.start_.begin(), h_ases.start_.end(),
                               indices.begin());

    // gather the remaining fields
    thrust::gather(indices.begin(), indices.end(), h_ases.end_.begin(),
                   h_ases_out.end_.begin());
    thrust::gather(indices.begin(), indices.end(), h_ases.strand.begin(),
                   h_ases_out.strand.begin());
    thrust::gather(indices.begin(), indices.end(), h_ases.core.begin(),
                   h_ases_out.core.begin());
    // swap
    h_ases.end_.swap(h_ases_out.end_);
    h_ases.strand.swap(h_ases_out.strand);
    h_ases.core.swap(h_ases_out.core);
}