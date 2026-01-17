#ifndef FASTQ_READER_H
#define FASTQ_READER_H
#include <string>
#include <htslib/kstring.h>
#include <htslib/hts.h>
#include <iostream>

class FastqPairReader {
public:
    struct Record {
        std::string header;
        std::string sequence;
        std::string quality;
        std::string plus;
    };

    struct FastqPair {
        Record r1;
        Record r2;
    };

    enum class ReadStatus {
        OK,
        END_OF_FILE,
        READ_ERROR,
    };
    FastqPairReader(const std::string& r1_path, const std::string& r2_path);
    ~FastqPairReader();
    ReadStatus next_record(FastqPair& pair);

private:
    std::string _r1_path_;
    std::string _r2_path_;
    htsFile *panel_r1 = nullptr;
    htsFile *panel_r2 = nullptr;
    // kstring_t has 3 fields: (1) l (length of string), m (allocated size for buffer), *s (string data).
    kstring_t line_r1 = KS_INITIALIZE;
    kstring_t line_r2 = KS_INITIALIZE;
    ReadStatus read_single_record(htsFile *fp, kstring_t &line, Record &rec);
    static std::string core_header(const std::string& header);
};

std::ostream& operator<<(std::ostream& os, const FastqPairReader::Record& rec);

std::ostream &operator<<(std::ostream &os, const FastqPairReader::FastqPair& pair);

#endif