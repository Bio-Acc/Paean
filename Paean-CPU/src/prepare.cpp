#include "gff.h"
#include "util.h"

#include <set>
#include <map>
#include <string>
#include <vector>
#include <cassert>
#include <utility>
#include <fstream>
#include <sstream>
#include <iostream>

#include <stdio.h>
#include <getopt.h>

void WriteGencodeAnnotation(std::string &gff_file, std::string &output_dir) {
    std::ifstream in_file(gff_file);
    if (!in_file.is_open()) {
        std::cerr << "Could not open this file: " << gff_file
                  << std::endl;
        exit(EXIT_FAILURE);
    }

    std::string ga_file = pathJoin(output_dir, "gene.annotation.gff3");
    std::ofstream out_file(ga_file);

    std::string line;
    while (std::getline(in_file, line)) {
        if (line.find("gene\t") != std::string::npos) {
            out_file << line << "\n";
        }
    }
    in_file.close();
    out_file.close();
}

uint32_t _getGeneLength(GffObj *gene) {
    uint32_t max_length = 0;

    for (int i = 0; i < gene->children.Count(); ++i) {
        GffObj *iso = gene->children[i];
        if (iso->isTranscript()) {
            uint32_t iso_length = 0;
            for (int j = 0; j < iso->exons.Count(); ++j) {
                GffExon *exon = iso->exons[j];
                iso_length += (exon->end - exon->start);
            }
            max_length = std::max(max_length, iso_length);
        }
    }
    return max_length;
}

void WriteLengthTable(std::string &gff_file, std::string &output_dir) {
    FILE *fileptr = fopen(gff_file.c_str(), "rb");
    if (!fileptr) {
        std::cerr << "Could not open this file: " << gff_file
                  << std::endl;
        exit(EXIT_FAILURE);
    }
    GffReader reader(fileptr);
    reader.readAll(true);

    std::string lt_file = pathJoin(output_dir, "length_table.csv");
    std::ofstream out_file(lt_file);

    if (!out_file.is_open()) {
        printf("error\n");
    }

    size_t nfeat = reader.gflst.Count();

    for (size_t i = 0; i < nfeat; ++i) {
        GffObj *f = reader.gflst[i];
        if (f->isGene()) {
            // gene id
            std::string gene_id = std::string(f->getGeneID());
            // gene length
            uint32_t gene_length = _getGeneLength(f);
            out_file << gene_id << "," << gene_length << "\n";
        }
    }
    out_file.close();
}

static void print_usage(FILE *fp) {
    fprintf(fp,
            "Usage: prepare [options...]\n"
            "Options:\n"
            "  -g,--gene            FILE     Gene annotation file with GFF3 format, required\n"
            "  -o,--output          FILE     Which directory you want to write the results to (default: ./)\n"
            "  -h,--help                     Print help information\n");
    exit(EXIT_SUCCESS);
}

int main(int argc, char **argv) {
    const char *const short_opts = "g:o:h";
    const option long_opts[] = {
        {"gff3", required_argument, nullptr, 'g'},
        {"output", required_argument, nullptr, 'o'},
        {"help", no_argument, nullptr, 'h'},
        {nullptr, no_argument, nullptr, 0}};

    std::string gff3_file, output_dir;

    int c;
    while ((c = getopt_long(argc, argv, short_opts, long_opts, nullptr)) >= 0) {
        switch (c) {
        case 'g':
            gff3_file = (char *)optarg;
            break;
        case 'o':
            output_dir = (char *)optarg;
            break;
        case 'h':
            print_usage(stdout);
            break;
        default:
            print_usage(stdout);
            break;
        }
    }

    if (gff3_file.empty()) {
        fprintf(stderr, "[prepare] require gene annotation file!\n");
        exit(EXIT_FAILURE);
    }

    if (!output_dir.empty()) {
        // if not exist, create
        createDirIfNotExists(output_dir.c_str());
    } else {
        // if no output_dir, use current dir
        output_dir = getPwd();
    }

    WriteGencodeAnnotation(gff3_file, output_dir);
    WriteLengthTable(gff3_file, output_dir);
}
