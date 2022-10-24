#include "parse.h"
#include "gff.h"
#include "util.h"
#include "fmt/core.h"
#include "htslib/sam.h"
#include "robin_hood.h"

#include <set>
#include <map>
#include <string>
#include <vector>
#include <cassert>
#include <utility>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>

// log file
extern std::ofstream log_file;

// global hash maps
UMAP ase_gid_map;
UMAP bin_gid_map;
UMAP_INT bin_len_map;

// temporary structure(AOS) for Read
struct t_Read {
    uint64_t start_ = 0;
    uint64_t end_ = 0;
    uint32_t length = 0;
    // uint8_t  strand;
    read_core_t core;
    bool is_null;

    t_Read() {
        is_null = false;
    }

    t_Read(bool is_null_) {
        is_null = is_null_;
    }
};
t_Read t_null_read(true);

// robin_hood_map
typedef robin_hood::unordered_map<size_t, 
                std::pair<t_Read, t_Read>> PAIR;
typedef robin_hood::unordered_map<size_t, t_Read> SINGLE;

// hash map of reference name's offset
robin_hood::unordered_map<std::string, uint32_t> offset_map;

inline bool is_null_read(t_Read &read) 
{ 
    return read.is_null;
}

uint64_t inline _offset(std::string &chr)
{
    auto it = offset_map.find(chr);
    if (it != offset_map.end()) {
        return it->second * refLength;
    }
    return invalidLength;
}

void LoadBinFromGff(h_Bins &h_bins, char *gff_file, char *csv_file)
{
    // load bin length table (length_table.csv)
    std::ifstream file(csv_file);
    if (file.fail()) {
        std::cerr << "Could not open this file: " << csv_file
                  << std::endl;
        exit(EXIT_FAILURE);
    }

    std::string header, bin_name, len_str;
    std::hash<std::string> hash_str;
    uint32_t len;
    
    while (std::getline(file, bin_name, ',')) {
        std::getline(file, len_str);
        len = std::stoi(len_str);
        bin_len_map.insert({hash_str(bin_name), len});
    }
    file.close();

    // load genes
    FILE *fileptr = fopen(gff_file, "rb");
    if (!fileptr) {
        std::cerr << "Could not open this file: " << gff_file
                  << std::endl;
        exit(EXIT_FAILURE);
    }
    GffReader reader(fileptr);
    reader.readAll(true);

    std::string gene_name, gene_id, chr_id, 
                strand, attrs;
    size_t nfeat = reader.gflst.Count();

    for (size_t i = 0; i < nfeat; ++i) {
        GffObj *f = reader.gflst[i];
        if (f->isGene()) {
            // chromosome id
            chr_id = std::string(f->getGSeqName());
            uint64_t offset = _offset(chr_id);
            // check offset
            if (offset == invalidLength) 
                continue;
                
            // gene name
            gene_name = std::string(f->getGeneName());
            gene_id = std::string(f->getGeneID());
            strand = f->strand;
            size_t hash_gid_t = hash_str(gene_id);
            attrs = join({gene_id, gene_name, chr_id, strand});
            bin_gid_map.emplace(hash_gid_t, attrs);
            
            // length
            size_t hash_name_t = hash_str(gene_name);
            auto it = bin_len_map.find(hash_name_t);
            if (it == bin_len_map.end())
                continue;
            uint32_t length = it->second;

            bin_core_t bin_core(hash_gid_t, length);
            h_bins.start_.emplace_back(f->start + offset);
            h_bins.end_.emplace_back(f->end + offset);
            h_bins.strand.emplace_back((f->strand == '+'));
            h_bins.core.emplace_back(bin_core);
// #define DEBUG
#ifdef DEBUG
            std::cout << "hash: " << hash_t << std::endl;
            std::cout << "gene hash: " << bin_core.gid_h << std::endl;
            std::cout << "gene start: " << f->start << std::endl;
            std::cout << "gene end: " << f->end << std::endl;
            std::cout << "gene strand: " << (f->strand == '+') << std::endl;
            std::cout << "gene length: " << length << std::endl;
#endif
        }
    }
}

// for this version, we unify the format (.csv) for different types of ase
void LoadAseFromCsv(h_ASEs &h_ases, char *csv_file)
{
    std::ifstream file(csv_file);
    if (file.fail()) {
        std::cerr << "Could not open this file: " << csv_file
                  << std::endl;
        exit(EXIT_FAILURE);
    }

    std::string ase_id, chr_id, strand, 
                out_donor, out_acceptor, in_junction;

    std::hash<std::string> hash_str;

    // skip header
    std::string header;
    std::getline(file, header);

    while (std::getline(file, ase_id, ',')) {
        std::getline(file, chr_id, ',');
        std::getline(file, strand, ',');
        std::getline(file, out_donor, ',');
        std::getline(file, out_acceptor, ',');
        std::getline(file, in_junction);
 
        // chromosome id
        uint64_t offset = _offset(chr_id);
        // check offset
        if (offset == invalidLength) 
            continue;
        
        // ase id
        size_t hash_t = hash_str(ase_id);
        // if duplicated, drop it
	    auto it = ase_gid_map.find(hash_t);
        if (it != ase_gid_map.end()) 
           continue;
        ase_gid_map.emplace(hash_t, ase_id);

        // strand
        uint8_t strand_ = strand.compare("+") == 0 ? 1 : 0;
        h_ases.strand.emplace_back(strand_);

        // coordinates
        ase_core_t ase_core(hash_t);
        int coordinateId = 0;

        // strip
        out_donor = out_donor.substr(0, out_donor.size() - 1);
        out_acceptor = out_acceptor.substr(0, out_acceptor.size() - 1);
        
        // out junction
        uint64_t donor = std::atoi(out_donor.c_str()) + offset;
        uint64_t acceptor = std::atoi(out_acceptor.c_str()) + offset;
        if (strand_) {
            ase_core.coordinates[coordinateId++] = donor;       
            ase_core.coordinates[coordinateId++] = acceptor;
        } else {
            ase_core.coordinates[coordinateId++] = acceptor;       
            ase_core.coordinates[coordinateId++] = donor;
        }
        
        // in junctions
        std::vector<uint64_t> coordinates;

        char *save_ptr = nullptr;
        char *split_str = strtok_r(const_cast<char*>(in_junction.c_str()), 
                                   "^-", &save_ptr);

        while (split_str) {
            if (isDigitstr(split_str)) {
                uint64_t coord = (uint64_t)strtol(split_str, nullptr, 10);
                coordinates.emplace_back(coord + offset); 
            }
            split_str = strtok_r(nullptr, "^-", &save_ptr);
        }
        if (!strand_) {
            std::reverse(coordinates.begin(), coordinates.end());
        }
        std::copy(coordinates.begin(), coordinates.end(),
                  ase_core.coordinates + coordinateId);
        ase_core.coordinateCountOut = coordinateId;
        ase_core.coordinateCountIn = coordinates.size();

        // ensure coordinateCountIn is an even number
        assert((ase_core.coordinateCountIn & 1) == 0);

        h_ases.core.emplace_back(ase_core);
        // in fact, we don't need to use start and end from now on.
        uint32_t coordinateCount = ase_core.coordinateCountOut + \
                                   ase_core.coordinateCountIn;
        uint64_t start_ = *std::min_element(ase_core.coordinates, 
                            ase_core.coordinates + coordinateCount);
        uint64_t end_ = *std::max_element(ase_core.coordinates, 
                            ase_core.coordinates + coordinateCount);
        h_ases.start_.emplace_back(start_);
        h_ases.end_.emplace_back(end_);

// #define DEBUG
#ifdef DEBUG
        uint32_t coordinateCount_ = ase_core.coordinateCountOut + \
                                    ase_core.coordinateCountIn;
        std::cout << "ase junction count: " << coordinateCount_
                  << std::endl;
        std::cout << "ase junctions: ";
        for (uint32_t i = 0; i < coordinateCount_; i++) {
            std::cout << ase_core.coordinates[i] << " ";
        }
        std::cout << std::endl;
#endif
    }
}

// -------------------------- pair-end --------------------------
// for read1 or read2 in a mate
static void _AddPairRead(PAIR &pairs, bam1_t *b, size_t hash_t, 
                         uint64_t offset, bool is_read1)
{
    /* we set both to false when there is only 
     * one read in a mate is chosen because there 
     * is no need to check if map `pairs` already 
     * contain hash_t, which saves time.
     */  

    read_core_t read_core;

    // extract junctions
    uint32_t *cigars = bam_get_cigar(b);
    uint32_t j;
    // skip flag `S` until `M`
    for (j = 0; j < b->core.n_cigar; j++) {
        if (bam_cigar_op(cigars[j]) == BAM_CMATCH)
            break;
    }
    // start acquiring junctions
    uint32_t prev = bam_cigar_oplen(cigars[j]);
    for (uint32_t i = j + 1; i < b->core.n_cigar; i++) {
        // flag == N
        uint32_t len = bam_cigar_oplen(cigars[i]);
        if (bam_cigar_op(cigars[i]) == BAM_CREF_SKIP) {
            read_core.junctions[read_core.junctionCount].start_ = prev;
            read_core.junctions[read_core.junctionCount].end_ = prev + len;
            read_core.junctionCount++;
        }
        prev = prev + len;
    }

    t_Read t_read;
    t_read.start_ = b->core.pos + offset + 1;
    t_read.end_ = bam_endpos(b) + offset + 1;
    t_read.length = b->core.l_qseq;
    // t_read.strand = (!bam_is_rev(b));
    t_read.core = read_core;

    if (is_read1) {
        // put read1 left
        auto it = pairs.find(hash_t);
        if (it != pairs.end()) {
            if (!is_null_read(it->second.first)) {
                if (read_core.junctionCount > 0) {
                    auto o_read = it->second.first;
                    std::copy(o_read.core.junctions,
                              o_read.core.junctions + o_read.core.junctionCount,
                              t_read.core.junctions + t_read.core.junctionCount);
                    t_read.core.junctionCount += o_read.core.junctionCount;
                    it->second.first = t_read;
                }
            } else {
                it->second.first = t_read;
            }
        } else {
            pairs.emplace(hash_t, std::make_pair(t_read, t_null_read));
        }
    } else {
        // put read2 right
        auto it = pairs.find(hash_t);
        if (it != pairs.end()) {
            if (!is_null_read(it->second.second)) {
                if (read_core.junctionCount > 0) {
                    auto o_read = it->second.second;
                    std::copy(o_read.core.junctions,
                              o_read.core.junctions + o_read.core.junctionCount,
                              t_read.core.junctions + t_read.core.junctionCount);
                    t_read.core.junctionCount += o_read.core.junctionCount;
                    it->second.second = t_read;
                }
            } else {
                it->second.second = t_read;
            }
        } else {
            pairs.emplace(hash_t, std::make_pair(t_null_read, t_read));
        }
    }
}

// merge read1 and read2 of a mate to one read
static void _MergePairRead(h_Reads &h_reads, h_Reads &h_nj_reads,
                           h_Junctions &h_junctions, h_Gaps &h_gaps,
                           PAIR &pairs, int max_gap)
{
    uint32_t numOfRead = 0, numOf_nj_Read = 0;

    for (auto &it : pairs) {
        t_Read &read1 = it.second.first;
        t_Read &read2 = it.second.second;

        // count number of reads with respect to with junctions
        // and without junctions
        if (read1.core.junctionCount > 0) {
            numOf_nj_Read++;
        } else {
            numOfRead++;
        }
        if (read2.core.junctionCount > 0) {
            numOf_nj_Read++;
        } else {
            numOfRead++;
        }

        // if one of these two reads is t_null_read
        // corresponds to single end
        if (is_null_read(read1) || is_null_read(read2)) {
            t_Read &read = is_null_read(read1) ? read2 : read1;
            if (read.core.junctionCount > 0) {
                h_reads.start_.emplace_back(read.start_);
                h_reads.end_.emplace_back(read.end_);
                // h_reads.strand.emplace_back(read.strand);
                // push junctions to junction table
                for (uint32_t j = 0; j < read.core.junctionCount; j++) {
                    h_junctions.start_.emplace_back(
                        read.start_ + read.core.junctions[j].start_ - 1);
                    h_junctions.end_.emplace_back(
                        read.start_ + read.core.junctions[j].end_);
                }
            } else {
                h_nj_reads.start_.emplace_back(read.start_);
                h_nj_reads.end_.emplace_back(read.end_);
                // h_nj_reads.strand.emplace_back(read.strand);
            }
            continue;
        }

        // merge start and end postions
        uint64_t start_ = std::min(read1.start_, read2.start_);
        uint64_t end_ = std::max(read1.end_, read2.end_);

        // if read1 and read2 have no overlap
        // here we consider the gap as a junction
        uint64_t gap_start = 0;
        uint64_t gap_end = 0;
        if (read1.end_ < read2.start_) {
            gap_start = read1.end_;
            gap_end = read2.start_;
        } else if (read2.end_ < read1.start_) {
            gap_start = read2.end_;
            gap_end = read1.start_;
        }
        // if inner distance is greater than `max_gap`, 
        // then discard this mate
        int dis = gap_end - gap_start;
        if (dis > 0) {
            if (dis > max_gap)
                continue;
            h_gaps.start_.emplace_back(gap_start);
            h_gaps.end_.emplace_back(gap_end);
        }

        // merge junctions
        // read1 == 0 and read2 == 0
        if (read1.core.junctionCount == 0 && 
            read2.core.junctionCount == 0) {
            h_nj_reads.start_.emplace_back(start_);
            h_nj_reads.end_.emplace_back(end_);
            // h_nj_reads.strand.emplace_back(read1.strand);
        } else {
            h_reads.start_.emplace_back(start_);
            h_reads.end_.emplace_back(end_);
            // h_reads.strand.emplace_back(read1.strand);
            
            read_core_t *read_core = nullptr;

            // read1 > 0 and read2 == 0
            if (read1.core.junctionCount > 0 && 
                read2.core.junctionCount == 0) {
                read_core = &read1.core;
            } 
            // read1 == 0 and read2 > 0
            else if (read1.core.junctionCount == 0 && 
                     read2.core.junctionCount > 0) {
                read_core = &read2.core;
            }
            // read1 > 0 and read2 > 0
            else if (read1.core.junctionCount > 0 && 
                     read2.core.junctionCount > 0) {
                // for union of two arrays
                std::set<std::pair<uint32_t, uint32_t>> set;

                // union of two junction arrays
                read_core_t read_core_merge;

                junction_t *junction1 = read1.core.junctions;
                junction_t *junction2 = read2.core.junctions;
                junction_t *junction = read_core_merge.junctions;

                // emplace read1 to set
                auto delta1 = uint32_t(read1.start_ - start_);
                for (uint32_t i1 = 0; i1 < read1.core.junctionCount; i1++) 
                    set.emplace(junction1[i1].start_+delta1, junction1[i1].end_+delta1);
                
                // emplace read2 to set
                auto delta2 = uint32_t(read2.start_ - start_);
                for (uint32_t i2 = 0; i2 < read2.core.junctionCount; i2++) 
                    set.emplace(junction2[i2].start_+delta2, junction2[i2].end_+delta2);

                // merge
                for (auto &pair : set) {
                    junction[read_core_merge.junctionCount].start_ = pair.first;
                    junction[read_core_merge.junctionCount].end_ = pair.second;
                    read_core_merge.junctionCount++;
                }
                read_core = &read_core_merge;
            }
            // h_reads.core.emplace_back(*read_core);
            // push junctions to junction table
            for (uint32_t j = 0; j < read_core->junctionCount; j++) {
                h_junctions.start_.emplace_back(
                    start_ + read_core->junctions[j].start_ - 1);
                h_junctions.end_.emplace_back(
                    start_ + read_core->junctions[j].end_);
            }
        }
    }

    log_file << "numOfRead\t" << numOfRead << "\n";
    log_file << "numOf_nj_Read\t" << numOf_nj_Read << "\n";
}
// --------------------------- end ------------------------------

// ------------------------- single-end -------------------------
static void _AddSingleRead(SINGLE &singles, bam1_t *b, 
                           size_t hash_t, uint64_t offset)
{
    read_core_t read_core;

    // extract junctions
    uint32_t *cigars = bam_get_cigar(b);
    uint32_t j;
    // skip flag `S` until `M`
    for (j = 0; j < b->core.n_cigar; j++) {
        if (bam_cigar_op(cigars[j]) == BAM_CMATCH)
            break;
    }
    // start acquiring junctions
    uint32_t prev = bam_cigar_oplen(cigars[j]);
    for (uint32_t i = j + 1; i < b->core.n_cigar; i++) {
        // flag == N
        uint32_t len = bam_cigar_oplen(cigars[i]);
        if (bam_cigar_op(cigars[i]) == BAM_CREF_SKIP) {
            read_core.junctions[read_core.junctionCount].start_ = prev;
            read_core.junctions[read_core.junctionCount].end_ = prev + len;
            read_core.junctionCount++;
        }
        prev = prev + len;
    }

    t_Read t_read;
    t_read.start_ = b->core.pos + offset + 1;
    t_read.end_ = bam_endpos(b) + offset + 1;
    t_read.length = b->core.l_qseq;
    // t_read.strand = (!bam_is_rev(b));
    t_read.core = read_core;

    auto it = singles.find(hash_t);
    if (it != singles.end()) {
        if (read_core.junctionCount > 0) {
            auto o_read = it->second;
            std::copy(o_read.core.junctions,
                      o_read.core.junctions + o_read.core.junctionCount,
                      t_read.core.junctions + t_read.core.junctionCount);
            t_read.core.junctionCount += o_read.core.junctionCount;
            it->second = t_read;
        }
    } else {
        singles.emplace(hash_t, t_read);
    }
}

// push read to vectors
static void _PushSingleRead(h_Reads &h_reads, h_Reads &h_nj_reads,
                            h_Junctions &h_junctions, SINGLE &singles)
{
    for (auto &it : singles) {
        t_Read &read = it.second;

        if (read.core.junctionCount > 0) {
            h_reads.start_.emplace_back(read.start_);
            h_reads.end_.emplace_back(read.end_);
            // h_reads.strand.emplace_back(read.strand);
            // push junctions to junction table
            for (uint32_t j = 0; j < read.core.junctionCount; j++) {
                h_junctions.start_.emplace_back(
                    read.start_ + read.core.junctions[j].start_ - 1);
                h_junctions.end_.emplace_back(
                    read.start_ + read.core.junctions[j].end_);
            }
        } else {
            h_nj_reads.start_.emplace_back(read.start_);
            h_nj_reads.end_.emplace_back(read.end_);
            // h_nj_reads.strand.emplace_back(read.strand);
        }
    }
}
// --------------------------- end ------------------------------

void LoadReadFromBam(h_Reads &h_reads, h_Reads &h_nj_reads,
                     h_Junctions &h_junctions, h_Gaps &h_gaps,
                     char *bam_file, int num_threads,
                     int max_gap, int mode)
{
    samFile *fp = sam_open(bam_file, "rb");
    if (!fp) {
        std::cerr << "Could not open this file: " << bam_file << std::endl;
        exit(EXIT_FAILURE);
    }
    
    // use multi-threading
    hts_set_threads(fp, num_threads);

    // header
    bam_hdr_t *header = sam_hdr_read(fp);
    if (!header) {
        std::cerr << "Could not read header for this file: " << bam_file
                  << std::endl;
        exit(EXIT_FAILURE);
    }

    // build hash map of reference name's offset
    for (uint32_t i = 0; i < (uint32_t)header->n_targets; i++) {
        auto name = std::string(header->target_name[i]);
        offset_map.emplace(name, i);
    }

    // parse
    uint64_t offset;
    
    // initialized structure of read
    bam1_t *b = bam_init1();

    /* From https://en.cppreference.com/w/cpp/utility/hash
     * For two different parameters k1 and k2 that 
     * are not equal, the probability that 
     * std::hash<Key>()(k1) == std::hash<Key>()(k2) 
     * should be very small, approaching 
     * 1.0/std::numeric_limits<std::size_t>::max()
     */
    std::hash<std::string> hash_str;

    auto is_invalid = [=](bam1_t *b) {
        bool invalid = ((b->core.tid < 0) ||
                        (b->core.flag & BAM_FUNMAP) ||
                        (b->core.flag & BAM_FSECONDARY));
        bool paired_end = ((mode == 2) ? true : false);
        bool is_paired = (b->core.flag & BAM_FPAIRED);
        return (invalid || (paired_end ^ is_paired));
    };

    // single-end
    if (mode == 1) {
        // temporary map for reads
        SINGLE singles;
        
        while (sam_read1(fp, header, b) >= 0) {
            // discard invalid reads
            if (is_invalid(b)) 
                continue;

            // we use tid as the offset of each read
            offset = b->core.tid * refLength;

            // construct unique name
            std::string read_name = bam_get_qname(b);
            std::string unique_name = fmt::format("{0};{1};{2}",
                            read_name, b->core.tid, b->core.pos);

            // we store integer instead of string
            // in order to save memory
            size_t hash_t = hash_str(unique_name);
            
            _AddSingleRead(singles, b, hash_t, offset);
        }

        // push read to vectors
        _PushSingleRead(h_reads, h_nj_reads, h_junctions, singles);
    }
    // pair-end
    else if (mode == 2) {
        // temporary map for fragments
        PAIR pairs;

        while (sam_read1(fp, header, b) >= 0) {
            // discard invalid reads
            if (is_invalid(b)) 
                continue;

            // we use tid as the offset of each read
            offset = b->core.tid * refLength;

            // hash qname
            std::string read_name = bam_get_qname(b);

            // remove suffix .1 and .2 or /1 and .2
            if (endsWith(read_name, ".1") || endsWith(read_name, ".2") ||
                endsWith(read_name, "/1") || endsWith(read_name, "/2")) {
                read_name = read_name.substr(0, read_name.size() - 2);
            }

            // because of multiple mapping, we must construct an unique name
            // construct unique read name
            std::string unique_name;
            if (b->core.isize > 0) {
                unique_name = fmt::format("{0};{1};{2};{3}",
                    read_name, b->core.tid, b->core.pos, b->core.mpos);
            } else {
                unique_name = fmt::format("{0};{1};{2};{3}",
                    read_name, b->core.tid, b->core.mpos, b->core.pos);
            }

            // we store integer instead of string
            // in order to save memory
            size_t hash_t = hash_str(unique_name);

            if (b->core.flag & BAM_FREAD1) {
                _AddPairRead(pairs, b, hash_t, offset, true);
            } else if (b->core.flag & BAM_FREAD2) {
                _AddPairRead(pairs, b, hash_t, offset, false);
            }
        }

        // merge pair-end reads
        _MergePairRead(h_reads, h_nj_reads, h_junctions, h_gaps,
                       pairs, max_gap);
    }

    assert(t_null_read.start_ == 0 && t_null_read.end_ == 0);
    // if no read with junctions
    if (!h_reads.start_.size()) {
        std::cerr << "This file has no valid reads with junctions: " << bam_file
                  << std::endl;
        exit(EXIT_FAILURE);
    }

// #define DEBUG
#ifdef DEBUG
    std::cout << "reads with junctions: " << std::endl;
    for (uint32_t itr = 0; itr < h_reads.start_.size(); itr++) {
        std::cout << "read start coordinate: " << h_reads.start_[itr] << std::endl;
        std::cout << "read end coordinate: " << h_reads.end_[itr] << std::endl;
        // std::cout << "read strand: " << (uint32_t)h_reads.strand[itr] << std::endl;
    }
    std::cout << "reads without junctions: " << std::endl;
    for (uint32_t itr = 0; itr < h_nj_reads.start_.size(); itr++) {
        std::cout << "read start coordinate: " << h_nj_reads.start_[itr] << std::endl;
        std::cout << "read end coordinate: " << h_nj_reads.end_[itr] << std::endl;
        // std::cout << "read strand: " << (uint32_t)h_nj_reads.strand[itr] << std::endl;
    }
    std::cout << "junctions: " << std::endl;
    std::cout << "junction count: " << h_junctions.start_.size() << std::endl;
    for (uint32_t itj = 0; itj < h_junctions.start_.size(); itj++) {
        std::cout << "junction start coordinate: " << h_junctions.start_[itj] << std::endl;
        std::cout << "junction end coordinate: " << h_junctions.end_[itj] << std::endl;
    }
    std::cout << "gaps: " << std::endl;
    std::cout << "gap count: " << h_gaps.start_.size() << std::endl;
    for (uint32_t itj = 0; itj < h_gaps.start_.size(); itj++) {
        std::cout << "gap start coordinate: " << h_gaps.start_[itj] << std::endl;
        std::cout << "gap end coordinate: " << h_gaps.end_[itj] << std::endl;
    }
#endif

    // free memory
    bam_destroy1(b);
    // close file
    sam_close(fp);
}
