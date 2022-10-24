#include "sort.cuh"

#include <thrust/scan.h>
#include <thrust/fill.h>
#include <thrust/replace.h>
#include <thrust/execution_policy.h>
#include <cub/device/device_scan.cuh>
#include <cub/device/device_reduce.cuh>
#include <cub/device/device_radix_sort.cuh>

void cubRadixSortKey(uint64_t *d_keys_in, uint64_t *d_keys_out,
                     uint32_t numOfEntry) {
    void *d_temp_storage = nullptr;
    size_t temp_storage_bytes = 0;

    cub::DeviceRadixSort::SortKeys(d_temp_storage, temp_storage_bytes,
                                   d_keys_in, d_keys_out, numOfEntry);
    // allocate temporary storage
    cudaMalloc(&d_temp_storage, temp_storage_bytes);
    // run sorting operation
    cub::DeviceRadixSort::SortKeys(d_temp_storage, temp_storage_bytes,
                                   d_keys_in, d_keys_out, numOfEntry);
}

void cubRadixSortInterval(d_Gaps &d_intervals_in, d_Gaps &d_intervals_out,
                          uint32_t numOfInterval) {
    // with junctions
    void *d_temp_storage = nullptr;
    size_t temp_storage_bytes = 0;
    // determine temporary device storage requirements
    cub::DeviceRadixSort::SortPairs(d_temp_storage, temp_storage_bytes,
                                    d_intervals_in.start_, d_intervals_out.start_,
                                    d_intervals_in.end_, d_intervals_out.end_,
                                    numOfInterval);
    // allocate temporary storage
    CUDA_SAFE_CALL(cudaMalloc(&d_temp_storage, temp_storage_bytes));
    // run sorting operation
    cub::DeviceRadixSort::SortPairs(d_temp_storage, temp_storage_bytes,
                                    d_intervals_in.start_, d_intervals_out.start_,
                                    d_intervals_in.end_, d_intervals_out.end_,
                                    numOfInterval);
    CUDA_SAFE_CALL(cudaDeviceSynchronize());
}

void cubRadixSortJunction(d_Junctions &d_junctions_in, d_Junctions &d_junctions_out,
                          uint32_t numOfJunction) {
    void *d_temp_storage = nullptr;
    size_t temp_storage_bytes = 0;
    // determine temporary device storage requirements
    cub::DeviceRadixSort::SortPairs(d_temp_storage, temp_storage_bytes,
                                    d_junctions_in.end_, d_junctions_out.end_,
                                    d_junctions_in.start_, d_junctions_out.start_,
                                    numOfJunction);
    // allocate temporary storage
    CUDA_SAFE_CALL(cudaMalloc(&d_temp_storage, temp_storage_bytes));
    // run sorting operation
    cub::DeviceRadixSort::SortPairs(d_temp_storage, temp_storage_bytes,
                                    d_junctions_in.end_, d_junctions_out.end_,
                                    d_junctions_in.start_, d_junctions_out.start_,
                                    numOfJunction);
    CUDA_SAFE_CALL(cudaDeviceSynchronize());

    d_temp_storage = nullptr;
    temp_storage_bytes = 0;
    // determine temporary device storage requirements
    cub::DeviceRadixSort::SortPairs(d_temp_storage, temp_storage_bytes,
                                    d_junctions_out.start_, d_junctions_in.start_,
                                    d_junctions_out.end_, d_junctions_in.end_,
                                    numOfJunction);
    // allocate temporary storage
    CUDA_SAFE_CALL(cudaMalloc(&d_temp_storage, temp_storage_bytes));
    // run sorting operation
    cub::DeviceRadixSort::SortPairs(d_temp_storage, temp_storage_bytes,
                                    d_junctions_out.start_, d_junctions_in.start_,
                                    d_junctions_out.end_, d_junctions_in.end_,
                                    numOfJunction);
    CUDA_SAFE_CALL(cudaDeviceSynchronize());

    CUDA_SAFE_CALL(cudaMemcpy(d_junctions_out.start_, d_junctions_in.start_,
                              sizeof(uint64_t) * numOfJunction,
                              cudaMemcpyDeviceToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(d_junctions_out.end_, d_junctions_in.end_,
                              sizeof(uint64_t) * numOfJunction,
                              cudaMemcpyDeviceToDevice));
}

d_Junctions thrustSegmentedScanJunction(d_Junctions &d_junctions_in, uint32_t &numOfJunction) {
    // segmented prefix sum
    uint32_t *d_counts_start, *d_counts_end, *d_counts;
    size_t jsize = sizeof(uint32_t) * numOfJunction;
    CUDA_SAFE_CALL(cudaMalloc((void **)&d_counts_start, jsize));
    CUDA_SAFE_CALL(cudaMalloc((void **)&d_counts_end, jsize));
    CUDA_SAFE_CALL(cudaMalloc((void **)&d_counts, jsize));

    // use thrust for pairwise segemented prefix sum
    thrust::fill(thrust::device, d_counts_start, d_counts_start + numOfJunction, 1);
    thrust::fill(thrust::device, d_counts_end, d_counts_end + numOfJunction, 1);
    thrust::inclusive_scan_by_key(thrust::device, d_junctions_in.start_,
                                  d_junctions_in.start_ + numOfJunction,
                                  d_counts_start, d_counts_start);
    thrust::inclusive_scan_by_key(thrust::device, d_junctions_in.end_,
                                  d_junctions_in.end_ + numOfJunction,
                                  d_counts_end, d_counts_end);
    thrust::transform(thrust::device, d_counts_start, d_counts_start + numOfJunction,
                      d_counts_end, d_counts, min_element<uint32_t>());

    // set flags
    uint32_t *d_flags;
    CUDA_SAFE_CALL(cudaMalloc((void **)&d_flags, jsize));
    thrust::replace_copy_if(thrust::device, d_counts, d_counts + numOfJunction,
                            d_flags, is_greater_than_one<uint32_t>(), 0);

    // compute offsets
    uint32_t *d_indices;
    CUDA_SAFE_CALL(cudaMalloc((void **)&d_indices, jsize));

    void *d_temp_storage = NULL;
    size_t temp_storage_bytes = 0;
    // determine temporary device storage requirements for inclusive prefix sum
    cub::DeviceScan::InclusiveSum(d_temp_storage, temp_storage_bytes, d_flags, d_indices,
                                  numOfJunction);
    // allocate temporary storage for inclusive prefix sum
    CUDA_SAFE_CALL(cudaMalloc(&d_temp_storage, temp_storage_bytes));
    // run inclusive prefix sum
    cub::DeviceScan::InclusiveSum(d_temp_storage, temp_storage_bytes, d_flags, d_indices,
                                  numOfJunction);
    CUDA_SAFE_CALL(cudaDeviceSynchronize());

    // calculate new numOfJunction
    uint32_t new_numOfJunction;
    CUDA_SAFE_CALL(cudaMemcpy(&new_numOfJunction, d_indices + numOfJunction - 1,
                              sizeof(uint32_t), cudaMemcpyDeviceToHost));

    d_Junctions d_junctions_out;
    CUDA_SAFE_CALL(
        cudaMalloc((void **)&d_junctions_out.start_, sizeof(uint64_t) * new_numOfJunction));
    CUDA_SAFE_CALL(
        cudaMalloc((void **)&d_junctions_out.end_, sizeof(uint64_t) * new_numOfJunction));
    CUDA_SAFE_CALL(
        cudaMalloc((void **)&d_junctions_out.count, sizeof(uint32_t) * new_numOfJunction));

    // compute number of junction block
    unsigned nJunctionBlock = DIV_UP(numOfJunction, blockSize);

    scatter_if<uint64_t><<<nJunctionBlock, blockSize>>>(
        d_junctions_in.start_, d_junctions_out.start_,
        d_indices, d_flags, numOfJunction);
    scatter_if<uint64_t><<<nJunctionBlock, blockSize>>>(
        d_junctions_in.end_, d_junctions_out.end_,
        d_indices, d_flags, numOfJunction);
    scatter_if<uint32_t><<<nJunctionBlock, blockSize>>>(
        d_counts, d_junctions_out.count,
        d_indices, d_flags, numOfJunction);
    CUDA_SAFE_CALL(cudaDeviceSynchronize());

    // update numOfJunction
    numOfJunction = new_numOfJunction;

    // free used memory
    CUDA_SAFE_CALL(cudaFree(d_counts_start));
    CUDA_SAFE_CALL(cudaFree(d_counts_end));
    CUDA_SAFE_CALL(cudaFree(d_counts));
    CUDA_SAFE_CALL(cudaFree(d_flags));
    CUDA_SAFE_CALL(cudaFree(d_indices));

    return d_junctions_out;
}

// sort bins by using cub library
void cubRadixSortBin(d_Bins &d_bins_in, d_Bins &d_bins_out,
                     h_Bins &h_bins, uint32_t numOfBin) {
    // indices on cpu
    auto *indices = new uint32_t[numOfBin];
    std::iota(indices, indices + numOfBin, 0);

    // indices on gpu
    uint32_t *d_indices_in;
    uint32_t *d_indices_out;
    CUDA_SAFE_CALL(
        cudaMalloc((void **)&d_indices_in, sizeof(uint32_t) * numOfBin));
    CUDA_SAFE_CALL(
        cudaMalloc((void **)&d_indices_out, sizeof(uint32_t) * numOfBin));
    CUDA_SAFE_CALL(
        cudaMemcpy(d_indices_in, indices, sizeof(uint32_t) * numOfBin,
                   cudaMemcpyHostToDevice));

    // with junctions
    void *d_temp_storage = nullptr;
    size_t temp_storage_bytes = 0;
    // determine temporary device storage requirements
    cub::DeviceRadixSort::SortPairs(d_temp_storage, temp_storage_bytes,
                                    d_bins_in.start_, d_bins_out.start_,
                                    d_indices_in, d_indices_out, numOfBin);
    // allocate temporary storage
    CUDA_SAFE_CALL(cudaMalloc(&d_temp_storage, temp_storage_bytes));
    // run sorting operation
    cub::DeviceRadixSort::SortPairs(d_temp_storage, temp_storage_bytes,
                                    d_bins_in.start_, d_bins_out.start_,
                                    d_indices_in, d_indices_out, numOfBin);
    CUDA_SAFE_CALL(cudaDeviceSynchronize());

    // compute number of thread block for bins
    unsigned nBinBlock = DIV_UP(numOfBin, blockSize);

    gather<uint64_t><<<nBinBlock, blockSize>>>(d_indices_out, d_bins_in.end_,
                                               d_bins_out.end_, numOfBin);
    gather<uint8_t><<<nBinBlock, blockSize>>>(d_indices_out, d_bins_in.strand,
                                              d_bins_out.strand, numOfBin);
    gather<bin_core_t><<<nBinBlock, blockSize>>>(d_indices_out, d_bins_in.core,
                                                 d_bins_out.core, numOfBin);
    CUDA_SAFE_CALL(cudaDeviceSynchronize());

    CUDA_SAFE_CALL(cudaMemcpy(h_bins.start_.data(), d_bins_out.start_,
                              numOfBin * sizeof(uint64_t), cudaMemcpyDeviceToHost));
    CUDA_SAFE_CALL(cudaMemcpy(h_bins.end_.data(), d_bins_out.end_,
                              numOfBin * sizeof(uint64_t), cudaMemcpyDeviceToHost));
    CUDA_SAFE_CALL(cudaMemcpy(h_bins.strand.data(), d_bins_out.strand,
                              numOfBin * sizeof(uint8_t), cudaMemcpyDeviceToHost));
    CUDA_SAFE_CALL(cudaMemcpy(h_bins.core.data(), d_bins_out.core,
                              numOfBin * sizeof(bin_core_t), cudaMemcpyDeviceToHost));

    // free used memory
    delete[] indices;
    CUDA_SAFE_CALL(cudaFree(d_indices_in));
    CUDA_SAFE_CALL(cudaFree(d_indices_out));
}

// sort ases by using cub library
void cubRadixSortASE(d_ASEs &d_ases_in, d_ASEs &d_ases_out,
                     uint32_t numOfASE) {
    // indices on cpu
    auto *indices = new uint32_t[numOfASE];
    std::iota(indices, indices + numOfASE, 0);

    // indices on gpu
    uint32_t *d_indices_in;
    uint32_t *d_indices_out;
    CUDA_SAFE_CALL(
        cudaMalloc((void **)&d_indices_in, sizeof(uint32_t) * numOfASE));
    CUDA_SAFE_CALL(
        cudaMalloc((void **)&d_indices_out, sizeof(uint32_t) * numOfASE));
    CUDA_SAFE_CALL(
        cudaMemcpy(d_indices_in, indices, sizeof(uint32_t) * numOfASE,
                   cudaMemcpyHostToDevice));

    // with junctions
    void *d_temp_storage = nullptr;
    size_t temp_storage_bytes = 0;
    // determine temporary device storage requirements
    cub::DeviceRadixSort::SortPairs(d_temp_storage, temp_storage_bytes,
                                    d_ases_in.start_, d_ases_out.start_,
                                    d_indices_in, d_indices_out, numOfASE);
    // allocate temporary storage
    CUDA_SAFE_CALL(cudaMalloc(&d_temp_storage, temp_storage_bytes));
    // run sorting operation
    cub::DeviceRadixSort::SortPairs(d_temp_storage, temp_storage_bytes,
                                    d_ases_in.start_, d_ases_out.start_,
                                    d_indices_in, d_indices_out, numOfASE);
    CUDA_SAFE_CALL(cudaDeviceSynchronize());

    // compute number of thread block for ases
    unsigned nASEBlock = DIV_UP(numOfASE, blockSize);

    gather<uint64_t><<<nASEBlock, blockSize>>>(d_indices_out, d_ases_in.end_,
                                               d_ases_out.end_, numOfASE);
    gather<uint8_t><<<nASEBlock, blockSize>>>(d_indices_out, d_ases_in.strand,
                                              d_ases_out.strand, numOfASE);
    gather<ase_core_t><<<nASEBlock, blockSize>>>(d_indices_out, d_ases_in.core,
                                                 d_ases_out.core, numOfASE);
    CUDA_SAFE_CALL(cudaDeviceSynchronize());

    // free used memory
    delete[] indices;
    CUDA_SAFE_CALL(cudaFree(d_indices_in));
    CUDA_SAFE_CALL(cudaFree(d_indices_out));
}

// cub reduce sum
void cubReduceSum(float *d_in, float *d_out, uint32_t num_items) {
    // determine temporary device storage requirements
    void *d_temp_storage = nullptr;
    size_t temp_storage_bytes = 0;
    int num_items_ = int(num_items);
    cub::DeviceReduce::Sum(d_temp_storage, temp_storage_bytes,
                           d_in, d_out, num_items_);
    // allocate temporary storage
    CUDA_SAFE_CALL(cudaMalloc(&d_temp_storage, temp_storage_bytes));
    // run sum-reduction
    cub::DeviceReduce::Sum(d_temp_storage, temp_storage_bytes,
                           d_in, d_out, num_items_);
    CUDA_SAFE_CALL(cudaDeviceSynchronize());
}
