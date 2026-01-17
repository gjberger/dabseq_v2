#include "fastq_reader.h"
#include <stdexcept>

using ReadStatus = FastqPairReader::ReadStatus;
using Record = FastqPairReader::Record;

FastqPairReader::FastqPairReader(const std::string &r1_path, const std::string &r2_path) {
    _r1_path_ = r1_path;
    _r2_path_ = r2_path;
    // Open R1/R2 in reading mode.
    // Note, HTSLib is written and C and requires
    // C strings, so we use `c_str()` to convert
    // std::string to a C string.
    panel_r1 = hts_open(r1_path.c_str(), "r");
    panel_r2 = hts_open(r2_path.c_str(), "r");

    // Throw error if fastq files don't open.
    if (!panel_r1 || !panel_r2) {
        // Clean up if only one is allocated, and the other isn't.
        if (panel_r1) hts_close(panel_r1);
        if (panel_r2) hts_close(panel_r2);
        throw std::runtime_error("Failed to open FASTQ files");
    }
}

FastqPairReader::~FastqPairReader() {
    ks_free(&line_r1);
    ks_free(&line_r2);
    // Free if not nullptr. HTSLib docs don't specify whether it's
    // safe to do on nullptr or not.
    if (panel_r1) hts_close(panel_r1); // Deallocate.
    if (panel_r2) hts_close(panel_r2); // Deallocate.
}

/*
    Extract the "core" Illumnia identifier from a FASTQ header,
    to validate R1/R2 read pairs.

    Illumina headers have the general form:
        @<instrument>:<run>:<flowcellID>:<lane>:<tile>:<x>:<y> <read>:<is_filtered>:<control>:<index>

    For example, an R1/R2 pair might look like:
        @LH00266:77:222WGNLT4:1:1101:47563:1028 1:N:0:GNAAGATC+AGTCGAAN
        @LH00266:77:222WGNLT4:1:1101:47563:1028 2:N:0:GNAAGATC+AGTCGAAN

    The chunk before the first space encodes metadata which is shared between R1/R2.
        @LH00266:77:222WGNLT4:1:1101:47563:1028

    The chunk after the first space is read-specific metadata:
         1:N:0:GNAAGATC+AGTCGAAN
         2:N:0:GNAAGATC+AGTCGAAN
        - "1" vs "2" → read number (R1 vs R2)
        - "N" vs "Y" → filter flag
        - control number → typically 0 for normal reads
        - index string → single or dual sample index (e.g. i7 + i5)
    
    This helper returns the substring from the beginning of the header up to
    (but not including) the space.

    Comparing the "core" of R1 and R2 lets us verify that the read is indeed a proper pair.
*/
std::string FastqPairReader::core_header(const std::string& header) {
    std::size_t space_pos = header.find(' ');

    // No space found: return whole header (eg @LH00266:77:222WGNLT4:4:1101:51131:1014)
    if (space_pos == std::string::npos) {
        return header;
    }

    // Return everything up to (but not including) the space.
    return header.substr(0, space_pos);
}

ReadStatus FastqPairReader::next_record(FastqPair &pair)
{
    ReadStatus read_one = read_single_record(panel_r1, line_r1, pair.r1);
    if (read_one == ReadStatus::END_OF_FILE) {
        return ReadStatus::END_OF_FILE;
    }
    
    if (read_one == ReadStatus::READ_ERROR) {
        return ReadStatus::READ_ERROR;
    }

    ReadStatus read_two = read_single_record(panel_r2, line_r2, pair.r2);
    // If R1 succeeded, but R2 did not, files are out of sync.
    if (read_two != ReadStatus::OK) {
        return ReadStatus::READ_ERROR;
    }

    std::string core_header_r1 = core_header(pair.r1.header);
    std::string core_header_r2 = core_header(pair.r2.header);
    if (core_header_r1 != core_header_r2) {
        return ReadStatus::READ_ERROR;
    }

    return ReadStatus::OK;
}

ReadStatus FastqPairReader::read_single_record(htsFile *file, kstring_t &line, Record &record) {
    // 1) HEADER
    int ret = hts_getline(file, '\n', &line);
    if (ret == -1) {
        return ReadStatus::END_OF_FILE;
    }

    if (ret <= -2) {
        return ReadStatus::READ_ERROR;
    }

    if (line.l == 0 || line.s[0] != '@') {
        return ReadStatus::READ_ERROR;
    }

    record.header.assign(line.s, line.l); // assign the k string to a std::string.

    // 2) SEQUENCE
    ret = hts_getline(file, '\n', &line);
    if (ret < 0) { // note, if END_OF_FILE occurs here, we have an incomplete read.
        return ReadStatus::READ_ERROR;
    }

    record.sequence.assign(line.s, line.l);

    // 3) PLUS LINE
    ret = hts_getline(file, '\n', &line);
    if (ret < 0) {
        return ReadStatus::READ_ERROR;
    }

    if (line.l == 0 || line.s[0] != '+'){
        return ReadStatus::READ_ERROR;
    }

    record.plus.assign(line.s, line.l);

    // 4) QUALITY
    ret = hts_getline(file, '\n', &line);

    if (ret < 0) {
        return ReadStatus::READ_ERROR;
    }

    record.quality.assign(line.s, line.l);

    // Sequence and quality lengths must match.
    if (record.sequence.size() != record.quality.size()) {
        return ReadStatus::READ_ERROR;
    }

    return ReadStatus::OK;
}

std::ostream &operator<<(std::ostream &os, const FastqPairReader::Record &rec) {
    os << rec.header << '\n'
       << rec.sequence << '\n'
       << rec.plus << '\n'
       << rec.quality;
    return os;
}

std::ostream &operator<<(std::ostream &os, const FastqPairReader::FastqPair &pair) {
    os << "R1\n" << pair.r1 << '\n'
       << "R2\n" << pair.r2;
    return os;
}