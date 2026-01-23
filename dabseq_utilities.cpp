#include "dabseq_utilities.h"
#include <fstream>
#include <sstream>
#include <stdexcept>

AntibodyPayloadResult extract_ab_payload_from_r2(const std::string &seq) {
    AntibodyPayloadResult result;
    result.valid = false;
    result.payload.clear();

    // Try Total-seq second option.
    std::size_t pos5 = find_with_mismatches(seq, H5_AB_HANDLE, 1); // tolerance of 1
    std::size_t pos3b = find_with_mismatches(seq, H3B_AB_HANDLE, 1);
    // MAYBE TRY CHECK FOR 3a?
    // DOESN'T NEED TO BE WHOLE H3B

    if (pos5 != std::string::npos && pos3b != std::string::npos) {
        std::size_t payload_start = pos5 + H5_AB_HANDLE.size();
        std::size_t payload_end = pos3b;

        // error check -> end is after start, and payload is large enough

        if (payload_end > payload_start) {
            result.payload = seq.substr(payload_start, payload_end - payload_start);
            result.valid = true;
            return result;
        }
    }

    std::size_t pos3a = find_with_mismatches(seq, H3A_AB_HANDLE, 1);
    // add length check
    if (pos3a != std::string::npos) {
        result.payload = seq.substr(0, pos3a);
        result.valid = true;
        return result;
    }

    return result; // neither pattern matched

}

std::size_t find_with_mismatches(const std::string &seq, const std::string &motif, int max_mismatches) {

    // If motif is empty or longer than the sequence, there is no match.
    if (motif.empty() || seq.size() < motif.size()) {
        return std::string::npos;
    }

    const std::size_t MOTIF_LENGTH = motif.size();
    const std::size_t SEQ_LENGTH = seq.size();

    for (std::size_t start = 0; start + MOTIF_LENGTH <= SEQ_LENGTH; start++) {
        int mismatches = 0;

        for (std::size_t i = 0; i < MOTIF_LENGTH; i++) {
            if (seq[start+i] != motif[i]) {
                mismatches++;
            }
            if (mismatches > max_mismatches) {
                break;
            }
        }

        if (mismatches <= max_mismatches) {
            return start;
        }
    }

    return std::string::npos;
}

ParsedBarcode parse_barcodes_from_r1(const FastqPairReader::Record &r1, const BarcodeIndex &barcodes)
{
    ParsedBarcode result = {"", "", false};

    const std::string &seq = r1.sequence;                 // efficient, read only reference to r1.sequence via seq

    std::size_t motif_pos = find_with_mismatches(seq, R1_START_MOTIF, 1); // hamming distance of 1
    // previously did .find, though this increases the number of barcodes found
    // flexibility in finding motif AND correcting barcodes makes more reads valid

    if (motif_pos == std::string::npos || motif_pos < 9)
        return result; // no motif or not enough bases before motif -> potentially don't need this, very pedantic

    // Raw observed barcodes from read.
    std::string barcode_first_half_observed = seq.substr(0, 9);
    std::string barcode_second_half_observed = seq.substr(motif_pos - 9, 9);

    // Canonical forms after Hamming correction.
    std::string barcode_first_half_corrected;
    std::string barcode_second_half_corrected;

    // Check if the barcodes can be mapped.
    bool barcode_first_half_mapped = barcodes.find_canonical_barcode(barcode_first_half_observed, barcode_first_half_corrected);
    bool barcode_second_half_mapped = barcodes.find_canonical_barcode(barcode_second_half_observed, barcode_second_half_corrected);

    if (!barcode_first_half_mapped || !barcode_second_half_mapped)
        return result; // at least one barcode isn't valid.

    result.valid = true;
    result.bc1 = std::move(barcode_first_half_corrected); // efficient transfer without copying
    result.bc2 = std::move(barcode_second_half_corrected);

    return result;
}

ParsedAntibody parse_antibody_from_r2(const FastqPairReader::Record& r2, const BarcodeIndex &antibody_barcodes) {
    ParsedAntibody result = {"", false};

    AntibodyPayloadResult pay = extract_ab_payload_from_r2(r2.sequence);
    if (!pay.valid || pay.payload.size() != 15) {
        return result; // not valid or not long enough
    }

    const std::string &barcode_observed = pay.payload;

    std::string barcode_corrected;

    bool barcode_mapped = antibody_barcodes.find_canonical_barcode(barcode_observed, barcode_corrected);

    if (!barcode_mapped) return result; // not in hamming dictionary.

    result.valid = true;
    result.barcode = std::move(barcode_corrected);

    return result;
}


std::unordered_map<std::string, std::string> load_antibody_name_map(const std::string& csv_path) {
    std::unordered_map<std::string, std::string> barcode_to_name;

    std::ifstream antibody_csv(csv_path);

    if (!antibody_csv) throw std::runtime_error("Failed to open antibody CSV: " + csv_path);
    
    std::string line;

    while(std::getline(antibody_csv, line)) {
        if(line.empty()) continue;

        std::size_t comma_pos = line.find(',');
        if (comma_pos == std::string::npos) throw std::runtime_error("Malformed antibody line (no comma)");

        std::string antibody_barcode = line.substr(0, comma_pos);
        std::string antibody_name = line.substr(comma_pos + 1);
        // ---- trim trailing CR / spaces / tabs from name ----
        while (!antibody_name.empty() &&
               (antibody_name.back() == '\r' ||
                antibody_name.back() == ' ' ||
                antibody_name.back() == '\t'))
        {
            antibody_name.pop_back();
        }

        // (optional) trim leading spaces/tabs from name
        std::size_t start = 0;
        while (start < antibody_name.size() &&
               (antibody_name[start] == ' ' || antibody_name[start] == '\t'))
        {
            ++start;
        }
        if (start > 0)
        {
            antibody_name.erase(0, start);
        }

        // (optional but nice) also trim trailing CR/space/tab on barcode
        while (!antibody_barcode.empty() &&
               (antibody_barcode.back() == '\r' ||
                antibody_barcode.back() == ' ' ||
                antibody_barcode.back() == '\t'))
        {
            antibody_barcode.pop_back();
        }

        barcode_to_name[antibody_barcode] = antibody_name;
    }

    return barcode_to_name;
}