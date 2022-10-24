#ifndef PAEAN_PARSE_H
#define PAEAN_PARSE_H

#include "bin.h"

// load bins from a GFF file
void LoadBinFromGff(h_Bins &, char*, char*);

// load ASEs from a csv file
void LoadAseFromCsv(h_ASEs &, char*);

// load reads from a bam file, including reads
// with junction and without junction.
// Moreover, we handle single-end and pair-end
// separately.
void LoadReadFromBam(h_Reads &, h_Reads &, 
                     h_Junctions &, h_Gaps &, char*,
                     int, int, int);

#endif  // PAEAN_PARSE_H