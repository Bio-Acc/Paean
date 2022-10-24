#pragma once

#include "bin.h"
#include "incgammabeta.cuh"

/* the kernel is used to count the tpm for each
 * bin. Note the d_bin_length is calculated by
 * parentheses algorithm.
 */
__global__ void gpu_count_tempTPM(d_Bins d_bins, uint32_t numOfBin,
                                  float *d_tempTPM) {
    uint32_t binId = blockDim.x * blockIdx.x + threadIdx.x;
    if (binId < numOfBin) {
        uint32_t length = d_bins.core[binId].length;
        if (length != 0) {
            d_tempTPM[binId] =
                float(d_bins.core[binId].readCount) / float(length);
        } else {
            d_tempTPM[binId] = 0;
        }
// #define DEBUG
#ifdef DEBUG
        printf("d_tempTPM: %f\n", d_tempTPM[binId]);
#endif
    }
}

/* the kernel is used to calculate an average tpm for
 * each bin. Note the d_tpmCount is calculated by
 * reduction algorithm.
 */
__global__ void gpu_count_TPM(d_Bins d_bins, uint32_t numOfBin,
                              float *d_tempTPM, float *d_tpmCount) {
    uint32_t binId = blockDim.x * blockIdx.x + threadIdx.x;

    if (binId < numOfBin) {
        // make sure d_tpmCounter is not zero
        if (*d_tpmCount == 0)
            return;
        // compute tpmCount for each gene
        float tpm = 1000000 * d_tempTPM[binId] / (*d_tpmCount);
        d_bins.core[binId].tpm = tpm;
// #define DEBUG
#ifdef DEBUG
        printf("d_tpmCount: %f\n", *d_tpmCount);
#endif
    }
}

/* the kernel is used to count the psi and confidence interval.
 */
__global__ void gpu_count_PSI(d_ASEs d_ases, uint32_t numOfASE,
                              ASEPsi *d_ase_psi) {
    uint32_t aseId = blockDim.x * blockIdx.x + threadIdx.x;
    float countIn, countOut, psi;
    float psi_ub, psi_lb, eps = 1.0e-5;

    if (aseId < numOfASE) {
        uint32_t nOut = d_ases.core[aseId].coordinateCountOut / 2;
        uint32_t nIn = d_ases.core[aseId].coordinateCountIn / 2;
        countIn = d_ase_psi[aseId].countIn / nIn;
        countOut = d_ase_psi[aseId].countOut / nOut;

        // compute pis
        if (countIn == 0 && countOut == 0) {
            psi = -1.0;
        } else {
            psi = (countIn) / (countIn + countOut);
        }

        // compute confidence interval
        psi_ub = 1 - invbetai(0.025, countOut, countIn + 1);
        psi_lb = 1 - invbetai(0.975, countOut + 1, countIn);

        if (fabs(countIn) < eps || fabs(countOut) < eps) {
            if (countIn + countOut >= 5) {
                psi_ub = 1;
                psi_lb = 1;
            } else {
                psi_ub = 1;
                psi_lb = 0;
            }
        }

        // store psi and confidence interval
        d_ase_psi[aseId].psi = psi;
        d_ase_psi[aseId].ciStart = psi_lb;
        d_ase_psi[aseId].ciEnd = psi_ub;

// #define DEBUG
#ifdef DEBUG
        if (aseId == 0) {
            for (int i = 0; i < numOfASE; i++) {
                printf("gid: %d, countIn: %.2f,  countOut: %.2f\n",
                       d_ase_psi[i].gid_h, d_ase_psi[i].countIn,
                       d_ase_psi[i].countOut);
            }
        }
#endif
    }
}