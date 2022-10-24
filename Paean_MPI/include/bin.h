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
#include <unordered_map>

// ase type
const uint32_t SE_Type   = 0;
const uint32_t A3SS_Type = 1;
const uint32_t A5SS_Type = 2;
const uint32_t RI_Type   = 3;
const uint32_t NSE_Type  = 4;

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

// for Assist, we use SOA
typedef struct {
    uint32_t *start_;
    uint32_t *end_;
} Assist;

// for read
struct read_core_t {
    // with junction
    uint32_t junctionCount = 0;
    junction_t junctions[junctionSize];
};

struct h_Reads {
    std::vector<uint64_t> start_;
    std::vector<uint64_t> end_;
    // std::vector<uint8_t> strand;
};

/* Note: there is no need to write constructor here
 * since `cudaMalloc` does not invoke constructor
 */
struct d_Reads {
    uint64_t *start_;
    uint64_t *end_;
    // uint8_t *strand;
};

/* here we need to build an extra junction table
 * for junctions of reads.
 */
struct h_Junctions {
    std::vector<uint64_t> start_;
    std::vector<uint64_t> end_;
};

struct d_Junctions {
    uint64_t *start_;
    uint64_t *end_;
    uint32_t *count;    // store how many times a junction occurs.
};

// only used for paired-end
struct h_Gaps {
    std::vector<uint64_t> start_;
    std::vector<uint64_t> end_;
};

struct d_Gaps {
    uint64_t *start_;
    uint64_t *end_;
};

// for bin
struct bin_core_t {
    size_t gid_h;
    uint32_t length;
    uint32_t readCount;
    float tpmCount;

    bin_core_t(size_t hash_, uint32_t length_) {
        readCount = tpmCount = 0;
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

struct d_Bins {
    uint64_t *start_;  //! absolute starting coordinate of bin
    uint64_t *end_;
    uint8_t *strand;

    bin_core_t *core;
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

struct d_ASEs {
    uint64_t *start_;
    uint64_t *end_;
    uint8_t *strand;

    ase_core_t *core;
};

typedef struct {
    size_t gid_h;
    bin_name_t bin_h;
    float countIn;
    float countOut;
    float psi;        // by sum(anchors)
    float ciStart;    // confidence interval
    float ciEnd;
} ASEPsi;

typedef std::unordered_map<size_t, std::string> UMAP;
typedef std::unordered_map<size_t, uint32_t> UMAP_INT;

struct pair_hash
{
    template <class T1, class T2>
    std::size_t operator() (const std::pair<T1, T2> &pair) const
    {
        return std::hash<T1>()(pair.first) ^ std::hash<T2>()(pair.second);
    }
};
typedef std::unordered_map<std::pair<size_t, size_t>, 
                           uint32_t, pair_hash> UMAP_PAIR_K;

// extern global hash maps
extern UMAP ase_gid_map;
extern UMAP bin_gid_map;
extern UMAP_INT bin_len_map;
extern UMAP_PAIR_K bin_read_map;

#endif // PAEAN_BIN_H