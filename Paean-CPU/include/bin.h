#ifndef PAEAN_BIN_H
#define PAEAN_BIN_H

#include <cmath>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <iomanip>
#include <vector>
#include <limits>
#include <string>

#include <robin_hood/robin_hood.h>

// ase type (int)
const uint32_t SE_Type = 0;
const uint32_t A3SS_Type = 1;
const uint32_t A5SS_Type = 2;
const uint32_t RI_Type = 3;
const uint32_t NSE_Type = 4;
const uint32_t numASETypes = 5;
// ase type (string)
const std::string aseTypes[] = {"SE", "A3SS", "A5SS", "RI", "NSE"};

// maximum field size
const int nameSize = 24;
const int gidSize = 96;
const int binNameSize = 5;
const int junctionSize = 10;
const int coordinateSize = 100;

const uint64_t refLength = 2ULL << 31;
const uint64_t invalidLength = refLength << 31;
const uint32_t refMask = (uint32_t)refLength - 1;

// kernel parameters
const int blockSize = 1024;

typedef struct {
    uint32_t start_ = 0;
    uint32_t end_ = 0;
} junction_t;

// for read
struct read_core_t {
    // with junction
    uint32_t junctionCount = 0;
    junction_t junctions[junctionSize];
};

struct h_Reads {
    std::vector<uint64_t> start_;
    std::vector<uint64_t> end_;
};

/* here we need to build an extra junction table
 * for junctions of reads.
 */
struct h_Junctions {
    std::vector<uint64_t> start_;
    std::vector<uint64_t> end_;
    std::vector<uint32_t> count;
};

// only used for paired-end
struct h_Gaps {
    std::vector<uint64_t> start_;
    std::vector<uint64_t> end_;
};

// for bin
struct bin_core_t {
    size_t gid_h;
    uint32_t length;
    uint32_t readCount;
    float tpm;

    bin_core_t(size_t hash_, uint32_t length_) {
        readCount = tpm = 0;
        gid_h = hash_;
        length = length_;
    }

    bin_core_t() {}
};

struct h_Bins {
    std::vector<uint64_t> start_;
    std::vector<uint64_t> end_;
    std::vector<uint8_t> strand;
    std::vector<bin_core_t> core;
};

// for ASE
struct bin_name_t {
    uint32_t binCount = 0;
    size_t bins[binNameSize];
};

struct ase_core_t {
    size_t gid_h;
    bin_name_t bin_h;
    uint32_t coordinateCountOut = 0;
    uint32_t coordinateCountIn = 0;
    uint64_t coordinates[coordinateSize];
    ase_core_t(size_t hash_) {
        gid_h = hash_;
        memset(coordinates, 0, sizeof(uint64_t) * coordinateSize);
    }

    ase_core_t() {}
};

struct h_ASEs {
    std::vector<uint64_t> start_;
    std::vector<uint64_t> end_;
    std::vector<uint8_t> strand;
    std::vector<ase_core_t> core;
};

typedef struct {
    size_t gid_h;
    bin_name_t bin_h;
    float countIn;
    float countOut;
    float psi;     // by sum(anchors)
    float ciStart; // confidence interval
    float ciEnd;
} ASEPsi;

struct h_Assist {
    std::vector<uint32_t> start_;
    std::vector<uint32_t> end_;
};

typedef robin_hood::unordered_map<size_t, std::string> UMAP;
typedef robin_hood::unordered_map<size_t, uint32_t> UMAP_INT;
typedef robin_hood::unordered_map<std::string, size_t> UMAP_REV;
// ase map
typedef robin_hood::unordered_map<size_t, h_ASEs> AMAP;

// extern global hash maps
extern UMAP ase_gid_map;
extern UMAP bin_gid_map;
extern UMAP_INT bin_len_map;

#endif // PAEAN_BIN_H