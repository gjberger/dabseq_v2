#ifndef DABSEQ_UTILITIES_H
#define DABSEQ_UTILITIES_H

#include <string>
#include "fastq_reader.h"
#include "barcode_index.h"

static const std::string R1_START_MOTIF = "GTACTCGCAGTAGTC";
static const std::string R1_END_MOTIF = "CTGTCTCTTATACACATCT";
static const std::string R2_END_MOTIF = "GACTACTGCGAGTAC";
static const std::string H5_AB_HANDLE = "TGACTACGCTACTCATGG";
static const std::string H3A_AB_HANDLE = "GCTTTAAGGCCGGTCCTAGC";
static const std::string H3B_AB_HANDLE = "GAGCCGATCTAGTATCTCAGTCG";

struct ParsedBarcode {
    std::string bc1;
    std::string bc2;
    bool valid;
};

struct AntibodyPayloadResult {
    std::string payload;
    bool valid;
};

struct ParsedAntibody {
    std::string barcode;
    bool valid;
};

ParsedAntibody parse_antibody_from_r2(const FastqPairReader::Record &r2, const BarcodeIndex &antibody_barcodes);

AntibodyPayloadResult extract_ab_payload_from_r2(const std::string &seq);

std::size_t find_with_mismatches(const std::string &seq, const std::string &motif, int max_mismatches);

ParsedBarcode parse_barcodes_from_r1(const FastqPairReader::Record &r1, const BarcodeIndex &barcodes);

std::unordered_map<std::string, std::string> load_antibody_name_map(const std::string &csv_path);

#endif // DABSEQ_UTILITIES_H