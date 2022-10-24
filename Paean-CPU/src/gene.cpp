#include <algorithm>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <chrono>
#include <ctime>

#include <omp.h>
#include <stdio.h>
#include <getopt.h>
#include <sys/stat.h>
#include <sys/sysinfo.h>

#include "bin.h"
#include "clock.h"
#include "parse.h"
#include "util.h"

#include "sort.h"
#include "bin_assign_kernel.h"
#include "bin_calc_kernel.h"

// log file
std::ofstream log_file;

void output_tpm(h_Bins &h_bins, uint32_t *h_readCounts,
                uint32_t *h_nj_readCounts,
                const std::string &output_dir,
                const std::string &file_prefix) {
    std::string bin_file = pathJoin(output_dir,
                                    file_prefix + "_gene.tsv");
    std::ofstream fout(bin_file);
    fout.setf(std::ios::fixed);

    // tpm
    fout << "gene_id\t"
         << "gene_name\t"
         << "reference\t"
         << "strand\t"
         << "start\t"
         << "end\t"
         << "reads_with_junction\t"
         << "reads_without_junction\t"
         << "total_reads\t"
         << "tpm\n";
    uint32_t size_ = h_bins.core.size();
    std::vector<std::string> attrs;
    for (uint32_t j = 0; j < size_; j++) {
        auto bin_start = h_bins.start_[j] & refMask;
        auto bin_end = h_bins.end_[j] & refMask;
        auto bin_core = h_bins.core[j];
        auto bin = bin_gid_map.find(bin_core.gid_h);
        assert(bin != bin_gid_map.end());
        attrs = split(bin->second);
        fout << attrs[0] << "\t"
             << attrs[1] << "\t"
             << attrs[2] << "\t"
             << attrs[3] << "\t"
             << bin_start << "\t"
             << bin_end << "\t"
             << h_readCounts[j] << "\t"
             << h_nj_readCounts[j] << "\t"
             << bin_core.readCount << "\t"
             << std::setprecision(3) << bin_core.tpm
             << "\n";
    }
    fout.close();
}

void output_psi(std::vector<ASEPsi> &h_ase_psi, uint32_t ase_type,
                const std::string &output_dir,
                const std::string &file_prefix) {
    std::string ase_file = pathJoin(output_dir,
                                    file_prefix + "_ase_" +
                                        aseTypes[ase_type] + ".tsv");
    std::ofstream fout(ase_file);
    fout.setf(std::ios::fixed);

    // print ASE type
    // fout << "ASE Type: " << aseTypes[ase_type] << "\n";

    // print header
    fout << "ase_id\t"
         << "gend_id\t"
         << "gene_name\t"
         << "count_in\t"
         << "count_out\t"
         << "psi\t"
         << "psi_lb\t"
         << "psi_ub\n";

    // uint32_t countIn, countOut;
    // countIn = countOut = 0;
    uint32_t numOfASE = h_ase_psi.size();

    for (uint32_t i = 0; i < numOfASE; i++) {
        // if (h_ase_psi[i].countIn) countIn++;
        // if (h_ase_psi[i].countOut) countOut++;

        auto gid = ase_gid_map.find(h_ase_psi[i].gid_h);
        fout << gid->second << "\t";

        std::stringstream ss_gene_id, ss_gene_name;
        uint32_t binCount = h_ase_psi[i].bin_h.binCount;
        if (binCount) {
            std::vector<std::string> attrs;
            for (uint32_t j = 0; j < binCount; j++) {
                auto bin = bin_gid_map.find(h_ase_psi[i].bin_h.bins[j]);
                if (bin == bin_gid_map.end()) {
                    ss_gene_id << "null";
                    ss_gene_name << "null";
                } else {
                    attrs = split(bin->second);
                    ss_gene_id << attrs[0];
                    ss_gene_name << attrs[1];
                }
                if (j != (binCount - 1)) {
                    ss_gene_id << ",";
                    ss_gene_name << ",";
                }
            }
        } else {
            ss_gene_id << "null";
            ss_gene_name << "null";
        }
        fout << ss_gene_id.str() << "\t";
        fout << ss_gene_name.str() << "\t";
        fout << std::setprecision(3) << h_ase_psi[i].countIn << "\t"
             << h_ase_psi[i].countOut << "\t"
             << h_ase_psi[i].psi << "\t"
             << h_ase_psi[i].ciStart << "\t"
             << h_ase_psi[i].ciEnd << "\n";
    }
    // fout << std::setprecision(3) << countIn << "\t"
    //      << std::setprecision(3) << countOut << "\n";
    fout.close();
}

void Run_tpm(h_Bins &h_bins, h_Reads h_reads, h_Reads h_nj_reads,
             h_Gaps &h_gaps, h_Gaps &h_nj_gaps,
             const std::string &output_dir,
             const std::string &file_prefix, bool paired_end) {
    uint32_t numOfBin = h_bins.start_.size();

    // sort gaps
    h_Gaps h_gap_pairs;
    // gaps
    h_gap_pairs.start_ = h_gaps.start_;
    h_gap_pairs.end_ = h_gaps.end_;
    thrust::sort(h_gaps.start_.begin(), h_gaps.start_.end());
    thrust::sort(h_gaps.end_.begin(), h_gaps.end_.end());
    thrust::sort_by_key(h_gap_pairs.start_.begin(),
                        h_gap_pairs.start_.end(),
                        h_gap_pairs.end_.begin());
    // nj_gaps
    h_Gaps h_nj_gap_pairs;
    h_nj_gap_pairs.start_ = h_nj_gaps.start_;
    h_nj_gap_pairs.end_ = h_nj_gaps.end_;
    thrust::sort(h_nj_gaps.start_.begin(), h_nj_gaps.start_.end());
    thrust::sort(h_nj_gaps.end_.begin(), h_nj_gaps.end_.end());
    thrust::sort_by_key(h_nj_gap_pairs.start_.begin(),
                        h_nj_gap_pairs.start_.end(),
                        h_nj_gap_pairs.end_.begin());

    // std::cout << "starting assign reads..." << std::endl;
    // calculate overlap
    std::vector<uint32_t> h_readCounts(numOfBin, 0);
    std::vector<uint32_t> h_nj_readCounts(numOfBin, 0);
    assign_read_kernel(h_bins, h_reads, h_gaps, h_gap_pairs,
                       h_readCounts.data(), paired_end);
    assign_read_kernel(h_bins, h_nj_reads, h_nj_gaps, h_nj_gap_pairs,
                       h_nj_readCounts.data(), paired_end);

    // std::cout << "starting count tpm..." << std::endl;

    std::vector<float> h_tempTPM(numOfBin, 0.0);
    count_tempTPM(h_bins, h_tempTPM.data());

    // calculate sum of tpms
    float h_tpmCount = thrust::reduce(h_tempTPM.begin(), h_tempTPM.end());

    // count tpm for each bin
    count_TPM(h_bins, h_tempTPM.data(), h_tpmCount);

    // output results to file
    output_tpm(h_bins, h_readCounts.data(), h_nj_readCounts.data(),
               output_dir, file_prefix);
}

void Run_psi(h_ASEs &h_ases, h_Bins &h_bins,
             h_Reads h_nj_reads, h_Junctions &h_junctions,
             uint32_t ase_type, const std::string &output_dir,
             const std::string &file_prefix) {
    uint32_t numOfASE = h_ases.start_.size();
    uint32_t numOfBin = h_bins.start_.size();

    // sort ases
    thrustSortASE(h_ases);

    // auxiliary array for ases
    h_Assist h_assist_ases;
    h_assist_ases.start_.resize(numOfBin, 0);
    h_assist_ases.end_.resize(numOfBin, 0);

    // assign ases to bins
    // std::cout << "starting assign ases..." << std::endl;
    assign_ASE_kernel(h_bins, h_ases, h_assist_ases);

    std::vector<ASEPsi> h_ase_psi(numOfASE);
    // assign reads to ases
    // std::cout << "starting assign reads to ases..." << std::endl;
    if (ase_type == RI_Type) {
        assign_read_RI_kernel(h_ases, h_nj_reads, h_junctions, h_ase_psi.data());
    } else {
        assign_read_ASE_kernel(h_ases, h_junctions, h_ase_psi.data());
    }

    // count psi and confidence interval
    // std::cout << "starting count psi and confidence interval..." << std::endl;
    count_PSI(h_ases, h_ase_psi.data());

    // output results to file
    output_psi(h_ase_psi, ase_type, output_dir, file_prefix);
}

void Run_paean(h_Bins &h_bins, h_Reads &h_reads, h_Reads &h_nj_reads,
               h_Gaps &h_gaps, h_Gaps &h_nj_gaps, h_Junctions &h_junctions,
               AMAP &h_ases_map, int gene_max_gap, int mode,
               const std::string &output_dir,
               const std::string &file_prefix) {
    if (!isdir(output_dir.c_str())) {
        std::cerr << "Could not find this directory: " << output_dir
                  << std::endl;
        exit(EXIT_FAILURE);
    }

    uint32_t numOfBin = h_bins.start_.size();
    uint32_t numOfRead = h_reads.start_.size();
    uint32_t numOf_nj_Read = h_nj_reads.start_.size();
    uint32_t numOfJunction = h_junctions.start_.size();

    log_file << "genes\t" << numOfBin << "\n";
    if (mode == 2) {
        log_file << "fragments with junctions\t" << numOfRead << "\n";
        log_file << "fragments without junctions\t" << numOf_nj_Read << "\n";
    } else {
        log_file << "reads with junctions\t" << numOfRead << "\n";
        log_file << "reads without junctions\t" << numOf_nj_Read << "\n";
    }
    log_file << "junctions\t" << numOfJunction << "\n";

    // std::cout << "starting sorting bins..." << std::endl;
    // sort bins
    thrustSortBin(h_bins);

    h_Reads h_reads_copy, h_nj_reads_copy;
    h_reads_copy.start_ = h_reads.start_;
    h_reads_copy.end_ = h_reads.end_;
    h_nj_reads_copy.start_ = h_nj_reads.start_;
    h_nj_reads_copy.end_ = h_nj_reads.end_;

    // std::cout << "starting sorting reads..." << std::endl;
    // sort read
    thrust::sort(h_reads.start_.begin(), h_reads.start_.end());
    thrust::sort(h_reads.end_.begin(), h_reads.end_.end());
    thrust::sort(h_nj_reads.start_.begin(), h_nj_reads.start_.end());
    thrust::sort(h_nj_reads.end_.begin(), h_nj_reads.end_.end());

    // std::cout << "starting sorting junctions..." << std::endl;
    // sort junctions
    thrustSortJunction(h_junctions);

    // std::cout << "starting counting junctions..." << std::endl;
    // count junctions
    h_Junctions h_cjunctions = thrustSegmentedScanJunction(h_junctions);

// #define DEBUG
#ifdef DEBUG
    numOfJunction = h_cjunctions.start_.size();
    printf("%u\n", numOfJunction);
    for (int i = 0; i < numOfJunction; i++)
        printf("%u %u %u %u\n", uint32_t(h_junctions.start_[i] / refLength),
               uint32_t(h_junctions.start_[i] & refMask),
               uint32_t(h_junctions.end_[i] & refMask),
               h_junctions.count_[i]);
#endif

    // std::cout << "starting computing tpm..." << std::endl;
    Run_tpm(h_bins, h_reads, h_nj_reads, h_gaps, h_nj_gaps,
            output_dir, file_prefix, (mode == 2));

    // std::cout << "starting computing psi..." << std::endl;
    // ASEs
    for (auto &it : h_ases_map) {
        Run_psi(it.second, h_bins, h_nj_reads, h_cjunctions,
                it.first, output_dir, file_prefix);
    }
}

#define DEAFULT_READ_MAX_GAP 500000
#define DEAFULT_GENE_MAX_GAP 500000

static void print_usage(FILE *fp) {
    const int num_threads = get_nprocs();

    fprintf(fp,
            "Usage: paean [options...]\n"
            "Options:\n"
            "  -g,--gene            FILE     Gene annotation file with GFF3 format, required\n"
            "  -l,--length          FILE     Gene length file with csv format, required\n"
            "  -x,--ase-types       STRING   Alternative splicing types, required\n"
            "  -y,--ase-files       FILE     ASE files with csv format, required\n"
            "  -r,--read            FILE     Read file with BAM format, required\n"
            "  -o,--output          FILE     Which directory you want to write the results to (default: ./)\n"
            "  -t,--thread          INT      Number of threads to use to parse BAM file (default: %d)\n"
            "  -a,--read-max-gap    INT      Max allowed gap for two reads of a mate (default: %d)\n"
            // "  -q,--gene-max-gap    INT      Max allowed gap for two genes for fusion detection (default: %d)\n"
            "  -m,--mode            INT      1 (single-end) or 2 (paired-end) mode\n"
            "  -s,--strandness      INT      Strand-specific read counting, 1 (stranded) or 2 (reversely stranded)\n"
            "  -M,--multi-mapping            Count Multi-mapping reads which are detected by 'NH' tag\n"
            "                                in the BAM/SAM file. All reported alignments will be counted\n"
            "                                for a multi-mapping read\n"
            "  -p,--primary                  Count primary alignments only\n"
            "  -h,--help                     Print help information\n",
            num_threads, DEAFULT_READ_MAX_GAP);
    exit(EXIT_SUCCESS);
}

int main(int argc, char **argv) {
    const char *const short_opts = "g:l:x:y:r:o:t:a:m:s:Mph";
    const option long_opts[] = {
        {"gene", required_argument, nullptr, 'g'},
        {"length", required_argument, nullptr, 'l'},
        {"ase-types", required_argument, nullptr, 'x'},
        {"ase-files", required_argument, nullptr, 'y'},
        {"read", required_argument, nullptr, 'r'},
        {"output", required_argument, nullptr, 'o'},
        {"thread", required_argument, nullptr, 't'},
        {"read-max-gap", required_argument, nullptr, 'a'},
        // {"gene-max-gap", required_argument, nullptr, 'q'},
        {"mode", required_argument, nullptr, 'm'},
        {"strandness", required_argument, nullptr, 's'},
        {"multi-mapping", no_argument, nullptr, 'M'},
        {"primary", no_argument, nullptr, 'p'},
        {"help", no_argument, nullptr, 'h'},
        {nullptr, no_argument, nullptr, 0}};

    struct GlobalArgs {
        char *bin_file;
        char *length_file;
        char *ase_types;
        char *ase_files;
        char *read_file;
        char *output_dir;
        int num_threads;
        int read_max_gap;
        int gene_max_gap;
        int mode;
        int strandness;
        bool multi_mapping;
        bool primary;

        GlobalArgs()
            : bin_file(nullptr), length_file(nullptr),
              ase_types(nullptr), ase_files(nullptr),
              read_file(nullptr), output_dir(nullptr),
              num_threads(get_nprocs()),
              read_max_gap(DEAFULT_READ_MAX_GAP),
              gene_max_gap(DEAFULT_GENE_MAX_GAP),
              mode(0), strandness(0), multi_mapping(false),
              primary(false) {}
    };
    GlobalArgs ga;

    int c;
    while ((c = getopt_long(argc, argv, short_opts, long_opts, nullptr)) >= 0) {
        switch (c) {
        case 'g':
            ga.bin_file = (char *)optarg;
            break;
        case 'l':
            ga.length_file = (char *)optarg;
            break;
        case 'x':
            ga.ase_types = (char *)optarg;
            break;
        case 'y':
            ga.ase_files = (char *)optarg;
            break;
        case 'r':
            ga.read_file = (char *)optarg;
            break;
        case 'o':
            ga.output_dir = (char *)optarg;
            break;
        case 't':
            ga.num_threads = std::stoi(optarg);
            break;
        case 'a':
            ga.read_max_gap = std::stoi(optarg);
            break;
        // case 'q':
        //     ga.gene_max_gap = std::stoi(optarg);
        //     break;
        case 'm':
            ga.mode = std::stoi(optarg);
            break;
        case 's':
            ga.strandness = std::stoi(optarg);
            break;
        case 'M':
            ga.multi_mapping = true;
            break;
        case 'p':
            ga.primary = true;
            break;
        case 'h':
            print_usage(stdout);
            break;
        default:
            print_usage(stdout);
            break;
        }
    }

    if (argc == 1) {
        print_usage(stdout);
    }

    // check input files
    if (!ga.bin_file || !ga.length_file || !ga.read_file) {
        fprintf(stderr, "[paean] require gene annotation file, length file, read file!\n");
        exit(EXIT_FAILURE);
    }

    // check ase
    if (!ga.ase_types || (ga.ase_types && !ga.ase_files)) {
        fprintf(stderr, "[paean] require the types and files of ASE to be specified!\n");
        exit(EXIT_FAILURE);
    }

    auto ase_types = split(std::string(ga.ase_types), ',');
    auto ase_files = split(std::string(ga.ase_files), ',');

    // check mathced size of ASE types and files
    if (ase_types.size() != ase_files.size()) {
        fprintf(stderr, "unmatched size of the types and files of ASE!\n");
        exit(EXIT_FAILURE);
    }

    // check ASE type name
    for (int i = 0; i < ase_types.size(); i++) {
        int idx = indexOf(aseTypes, numASETypes, ase_types[i]);
        if (idx == -1) {
            fprintf(stderr, "wrong ASE type!\n");
            exit(EXIT_FAILURE);
        }
    }

    // check mode
    if (ga.mode != 1 && ga.mode != 2) {
        fprintf(stderr, "please specify one mode, 1(single-end) or 2(pair-end)!\n");
        exit(EXIT_FAILURE);
    }

    // check strandness
    if (ga.strandness != 0 && ga.strandness != 1 && ga.strandness != 2) {
        fprintf(stderr,
                "please specify correct strandness, 1 (stranded) or"
                " 2 (reversely stranded)!\n");
        exit(EXIT_FAILURE);
    }

    std::string output_dir;
    if (ga.output_dir) {
        output_dir = std::string(ga.output_dir);
        // if not exist, create
        createDirIfNotExists(output_dir.c_str());
    } else {
        // if no output_dir, use current dir
        output_dir = getPwd();
    }

    // set number of threads
    omp_set_num_threads(ga.num_threads);

    // extract file prefix
    std::string file_prefix = extractFilePrefix(ga.read_file);

    // open log file
    std::string log_name = pathJoin(output_dir, file_prefix + "_log.txt");
    log_file.open(log_name);

    // bins and ases
    h_Bins h_bins;
    // h_ASEs h_se_ases, h_a3ss_ases, h_a5ss_ases;
    // reads
    h_Reads h_reads, h_nj_reads;
    h_Gaps h_gaps, h_nj_gaps;
    h_Junctions h_junctions;

    using namespace std::chrono;
    auto t0 = high_resolution_clock::now();

    // load reads
    // std::cout << "loading reads..." << std::endl;
    int num_jreads = 0;
    LoadReadFromBam(h_reads, h_nj_reads, h_gaps,
                    h_nj_gaps, h_junctions, ga.read_file,
                    ga.num_threads, ga.read_max_gap,
                    ga.mode, ga.strandness, ga.multi_mapping,
                    ga.primary);

    // load bins
    // std::cout << "loading bins..." << std::endl;
    LoadBinFromGff(h_bins, ga.bin_file, ga.length_file);

    // load ases
    AMAP h_ases_map;
    for (int i = 0; i < ase_types.size(); i++) {
        size_t idx = indexOf(aseTypes, numASETypes, ase_types[i]);
        h_ASEs h_ases;
        LoadAseFromCsv(h_ases, &ase_files[i][0]);
        h_ases_map.emplace(idx, h_ases);
    }

    auto t1 = high_resolution_clock::now();

    // handling paean
    // std::cout << "start kernel program..." << std::endl;

    Run_paean(h_bins, h_reads, h_nj_reads, h_gaps, h_nj_gaps,
              h_junctions, h_ases_map, ga.gene_max_gap,
              ga.mode, output_dir, file_prefix);

    auto t2 = high_resolution_clock::now();

    // set precision
    log_file.setf(std::ios::fixed);
    log_file.precision(1);

    log_file << "load time\t"
             << duration_cast<duration<double>>(t1 - t0).count()
             << "s\n";
    log_file << "computation time\t"
             << duration_cast<duration<double>>(t2 - t1).count()
             << "s\n";
    log_file.close();
}
