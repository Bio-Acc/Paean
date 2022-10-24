#include "parse.h"
#include "gff.h"
#include "util.h"
#include "sam.h"

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
    uint8_t strand;
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

// vector
typedef std::pair<t_Read, t_Read> FRAGMENT;

// hash map of reference name's offset
UMAP_REV offset_map;

inline bool is_null_read(t_Read &read) {
    return read.is_null;
}

inline uint64_t _offset(std::string &chr) {
    auto it = offset_map.find(chr);
    if (it != offset_map.end()) {
        return it->second * refLength;
    }
    return invalidLength;
}

void LoadBinFromGff(h_Bins &h_bins, char *gff_file, char *csv_file) {
    // load bin length table (length_table.csv)
    std::ifstream file(csv_file);
    if (file.fail()) {
        std::cerr << "Could not open this file: " << csv_file
                  << std::endl;
        exit(EXIT_FAILURE);
    }

    std::string header, bin_id, len_str;
    std::hash<std::string> hash_str;
    uint32_t len;

    while (std::getline(file, bin_id, ',')) {
        std::getline(file, len_str);
        len = std::stoi(len_str);
        bin_len_map.insert({hash_str(bin_id), len});
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
            gene_id = std::string(f->getGeneID());
            if (!f->getGeneName()) {
                gene_name = gene_id;
            } else {
                gene_name = std::string(f->getGeneName());
            }
            strand = f->strand;
            size_t hash_gid_t = hash_str(gene_id);
            attrs = join({gene_id, gene_name, chr_id, strand});
            bin_gid_map.emplace(hash_gid_t, attrs);

            // length
            auto it = bin_len_map.find(hash_gid_t);
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
void LoadAseFromCsv(h_ASEs &h_ases, char *csv_file) {
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
        char *split_str = strtok_r(const_cast<char *>(in_junction.c_str()),
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
        uint32_t coordinateCount = ase_core.coordinateCountOut +
                                   ase_core.coordinateCountIn;
        uint64_t start_ = *std::min_element(ase_core.coordinates,
                                            ase_core.coordinates + coordinateCount);
        uint64_t end_ = *std::max_element(ase_core.coordinates,
                                          ase_core.coordinates + coordinateCount);
        h_ases.start_.emplace_back(start_);
        h_ases.end_.emplace_back(end_);

// #define DEBUG
#ifdef DEBUG
        uint32_t coordinateCount_ = ase_core.coordinateCountOut +
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

static t_Read _ConvertToRead(bam1_t *b, uint64_t offset) {
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
    t_read.strand = (!bam_is_rev(b));
    t_read.core = read_core;

    return t_read;
}

// -------------------------- pair-end --------------------------

// merge read1 and read2 of a mate to one read
static void _MergePairRead(h_Reads &h_reads, h_Reads &h_nj_reads,
                           h_Gaps &h_gaps, h_Gaps &h_nj_gaps,
                           h_Junctions &h_junctions, FRAGMENT &frag,
                           int max_gap, int strandness) {
    t_Read &read1 = frag.first;
    t_Read &read2 = frag.second;

    // if one of these two reads is t_null_read
    // corresponds to single end
    if (is_null_read(read1) || is_null_read(read2)) {
        t_Read &read = is_null_read(read1) ? read2 : read1;
        if (read.core.junctionCount > 0) {
            h_reads.start_.emplace_back(read.start_);
            h_reads.end_.emplace_back(read.end_);
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
        }
        return;
    }

    // check strand
    if ((strandness == 1 && read1.strand != 1 && read2.strand != 0) ||
        (strandness == 2 && read1.strand != 0 && read2.strand != 1)) {
        return;
    }

    // merge start and end postions
    uint64_t start_ = std::min(read1.start_, read2.start_);
    uint64_t end_ = std::max(read1.end_, read2.end_);

    // if read1 and read2 have no overlap
    // here we consider the impact of gap
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
    if (dis > max_gap)
        return;

    // merge junctions
    // read1 == 0 and read2 == 0
    if (read1.core.junctionCount == 0 &&
        read2.core.junctionCount == 0) {
        h_nj_reads.start_.emplace_back(start_);
        h_nj_reads.end_.emplace_back(end_);

        if (dis > 0) {
            h_nj_gaps.start_.emplace_back(gap_start);
            h_nj_gaps.end_.emplace_back(gap_end);
        }
    } else {
        h_reads.start_.emplace_back(start_);
        h_reads.end_.emplace_back(end_);

        if (dis > 0) {
            h_gaps.start_.emplace_back(gap_start);
            h_gaps.end_.emplace_back(gap_end);
        }

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
                set.emplace(junction1[i1].start_ + delta1, junction1[i1].end_ + delta1);

            // emplace read2 to set
            auto delta2 = uint32_t(read2.start_ - start_);
            for (uint32_t i2 = 0; i2 < read2.core.junctionCount; i2++)
                set.emplace(junction2[i2].start_ + delta2, junction2[i2].end_ + delta2);

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
// --------------------------- end ------------------------------

// ------------------------- single-end -------------------------

// push read to vectors
static void _PushSingleRead(h_Reads &h_reads, h_Reads &h_nj_reads,
                            h_Junctions &h_junctions, t_Read &read) {
    if (read.core.junctionCount > 0) {
        h_reads.start_.emplace_back(read.start_);
        h_reads.end_.emplace_back(read.end_);
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
    }
}
// --------------------------- end ------------------------------

void LoadReadFromBam(h_Reads &h_reads, h_Reads &h_nj_reads,
                     h_Gaps &h_gaps, h_Gaps &h_nj_gaps,
                     h_Junctions &h_junctions, char *bam_file,
                     int num_threads, int max_gap, int mode,
                     int strandness, bool multi_mapping,
                     bool primary) {
    samFile *fp = sam_open(bam_file, "rb");
    if (!fp) {
        std::cerr << "Could not open the file: " << bam_file << std::endl;
        exit(EXIT_FAILURE);
    }

    // use multi-threading
    hts_set_threads(fp, num_threads);

    // header
    bam_hdr_t *header = sam_hdr_read(fp);
    if (!header) {
        std::cerr << "Could not read header for the file: " << bam_file
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

    // read the first sam line
    if (sam_read1(fp, header, b) < 0) {
        std::cerr << "Could not read alignments for the file: " << bam_file
                  << std::endl;
        exit(EXIT_FAILURE);
    }

    // check if 'NH' tag is in the BAM file
    uint8_t *tag = bam_aux_get(b, "NH");
    if (!multi_mapping && tag == NULL) {
        multi_mapping = true;
        primary = false;
    }

    bool paired_end = ((mode == 2) ? true : false);

    auto is_invalid = [=](bam1_t *b) {
        bool invalid = ((b->core.tid < 0) ||
                        (b->core.flag & BAM_FUNMAP));
        bool is_paired = (b->core.flag & BAM_FPAIRED);
        invalid = (invalid || (paired_end ^ is_paired));

        if (invalid)
            return invalid;

        if (multi_mapping) {
            if (primary) {
                invalid = (b->core.flag & BAM_FSECONDARY);
            } else {
                if ((b->core.flag & BAM_FSECONDARY) &&
                    !(b->core.flag & BAM_FPROPER_PAIR)) {
                    invalid = true;
                }
            }
        } else {
            if (b->core.flag & BAM_FSECONDARY) {
                invalid = true;
            } else {
                uint8_t *tag = bam_aux_get(b, "NH");
                int nh = bam_aux2i(tag);
                if (nh > 1) {
                    invalid = true;
                }
            }
        }

        return invalid;
    };

    if (paired_end) {
        do {
            // discard invalid reads
            if (is_invalid(b))
                continue;

            FRAGMENT frag;

            // we use tid as the offset of each read
            offset = b->core.tid * refLength;

            if (b->core.flag & BAM_FMUNMAP) {
                // if only one read mapped
                t_Read t_read = _ConvertToRead(b, offset);

                if (b->core.flag & BAM_FREAD1) {
                    frag.first = t_read;
                    frag.second = t_null_read;
                } else {
                    frag.first = t_null_read;
                    frag.second = t_read;
                }
            } else {
                // two reads mapped
                t_Read f_read = _ConvertToRead(b, offset);

                // extract the second read
                int status = sam_read1(fp, header, b);
                assert(status >= 0);

                uint64_t s_offset = b->core.tid * refLength;
                // filter paired reads mapped to different chrs
                if (s_offset != offset)
                    continue;

                t_Read s_read = _ConvertToRead(b, s_offset);

                if (b->core.flag & BAM_FREAD2) {
                    frag.first = f_read;
                    frag.second = s_read;
                } else {
                    frag.first = s_read;
                    frag.second = f_read;
                }
            }

            // merge pair-end reads
            _MergePairRead(h_reads, h_nj_reads, h_gaps, h_nj_gaps,
                           h_junctions, frag, max_gap, strandness);

        } while (sam_read1(fp, header, b) >= 0);
    } else {
        // single-end
        do {
            // discard invalid reads
            if (is_invalid(b))
                continue;

            // we use tid as the offset of each read
            offset = b->core.tid * refLength;

            t_Read t_read = _ConvertToRead(b, offset);

            // push read to vectors
            _PushSingleRead(h_reads, h_nj_reads, h_junctions, t_read);

        } while (sam_read1(fp, header, b) >= 0);
    }

    // if no read with junctions
    if (!h_reads.start_.size()) {
        std::cerr << "The file has no valid reads with junctions: " << bam_file
                  << std::endl;
        exit(EXIT_FAILURE);
    }

// #define DEBUG
#ifdef DEBUG
    std::cout << "reads with junctions: " << h_reads.start_.size() << std::endl;
    for (uint32_t itr = 0; itr < h_reads.start_.size(); itr++) {
        std::cout << "read start coordinate: " << h_reads.start_[itr] << std::endl;
        std::cout << "read end coordinate: " << h_reads.end_[itr] << std::endl;
    }
    std::cout << "reads without junctions: " << h_nj_reads.start_.size() << std::endl;
    for (uint32_t itr = 0; itr < h_nj_reads.start_.size(); itr++) {
        std::cout << "read start coordinate: " << h_nj_reads.start_[itr] << std::endl;
        std::cout << "read end coordinate: " << h_nj_reads.end_[itr] << std::endl;
    }
    std::cout << "junctions: " << h_junctions.start_.size() << std::endl;
    for (uint32_t itj = 0; itj < h_junctions.start_.size(); itj++) {
        std::cout << "junction start coordinate: " << h_junctions.start_[itj] << std::endl;
        std::cout << "junction end coordinate: " << h_junctions.end_[itj] << std::endl;
    }
    std::cout << "gaps: " << h_gaps.start_.size() << std::endl;
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
