#ifndef PAEAN_BIN_CALC_KERNEL_H
#define PAEAN_BIN_CALC_KERNEL_H

#include "bin.h"
#include "incgammabeta.h"

/* the kernel is used to count the tpm for each
 * bin. Note the d_bin_length is calculated by
 * parentheses algorithm.
 */
void count_tempTPM(h_Bins &h_bins, float *h_tempTPM) {
    uint32_t numOfBin = h_bins.start_.size();

#pragma omp parallel for
    for (uint32_t binId = 0; binId < numOfBin; binId++) {
        uint32_t length = h_bins.core[binId].length;
        if (length != 0) {
            h_tempTPM[binId] =
                float(h_bins.core[binId].readCount) / float(length);
        } else {
            h_tempTPM[binId] = 0;
        }
// #define DEBUG
#ifdef DEBUG
        printf("h_tempTPM: %f\n", h_tempTPM[binId]);
#endif
    }
}

/* the kernel is used to calculate an average tpm for
 * each bin. Note the h_tpmCount is calculated by
 * reduction algorithm.
 */
void count_TPM(h_Bins &h_bins, float *h_tempTPM,
               float h_tpmCount) {
    uint32_t numOfBin = h_bins.start_.size();

#pragma omp parallel for
    for (uint32_t binId = 0; binId < numOfBin; binId++) {
        if (h_tpmCount == 0)
            continue;
        // compute tpmCount for each gene
        float tpm = 1000000 * h_tempTPM[binId] / (h_tpmCount);
        h_bins.core[binId].tpm = tpm;
    }
}

/* the kernel is used to count the psi and confidence interval.
 */
void count_PSI(h_ASEs &h_ases, ASEPsi *h_ase_psi) {
    uint32_t numOfASE = h_ases.start_.size();

#pragma omp parallel for
    for (uint32_t aseId = 0; aseId < numOfASE; aseId++) {
        float countIn, countOut, psi;
        float psi_ub, psi_lb, eps = 1.0e-5;

        uint32_t nOut = h_ases.core[aseId].coordinateCountOut / 2;
        uint32_t nIn = h_ases.core[aseId].coordinateCountIn / 2;
        countIn = h_ase_psi[aseId].countIn / nIn;
        countOut = h_ase_psi[aseId].countOut / nOut;

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
        h_ase_psi[aseId].psi = psi;
        h_ase_psi[aseId].ciStart = psi_lb;
        h_ase_psi[aseId].ciEnd = psi_ub;

// #define DEBUG
#ifdef DEBUG
        for (uint32_t aseId = 0; aseId < numOfASE; aseId++) {
            printf("gid: %d, countIn: %.2f,  countOut: %.2f\n",
                   h_ase_psi[aseId].gid_h, h_ase_psi[aseId].countIn,
                   h_ase_psi[aseId].countOut);
        }
#endif
    }
}

#endif // PAEAN_BIN_CALC_KERNEL_H