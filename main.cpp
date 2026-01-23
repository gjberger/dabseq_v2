#include <iostream>
#include "fastq_reader.h"
#include "barcode_index.h"
#include "dabseq_utilities.h"
#include <fstream>
#include <filesystem>

namespace fs = std::filesystem;

int main(int argc, char **argv)
{

    if (argc != 5)
    {
        std::cerr << "Usage: " << argv[0] << " R1.fastq[.gz] R2.fastq[.gz] cell_barcodes.csv antibody_barcodes.csv\n";
        return 1;
    }

    const std::string r1_path = argv[1];
    const std::string r2_path = argv[2];
    const std::string cell_barcodes_csv = argv[3];
    const std::string antibody_barcodes_csv = argv[4];

    try
    {
        // Load cell barcode whitelist.
        BarcodeIndex cell_barcode_set(cell_barcodes_csv);
        // Load antibody barcode whitelist.
        BarcodeIndex antibody_barcode_set(antibody_barcodes_csv);

        std::unordered_map<std::string, std::string> antibody_barcode_to_name = load_antibody_name_map(antibody_barcodes_csv);

        FastqPairReader reader(r1_path, r2_path);
        FastqPairReader::FastqPair pair;

        int printed_pairs = 0;
        const int max_pairs = 10000;
        std::size_t num_with_barcodes = 0; // measure effectiveness of hamming dictionary
        std::size_t num_with_ab_payload = 0;
        std::size_t num_with_both = 0;

        //                  cell_id     antibody, cnt
        std::unordered_map<std::string, std::unordered_map<std::string, std::size_t>> counts;

        while (printed_pairs < max_pairs)
        {
            FastqPairReader::ReadStatus status = reader.next_record(pair);

            if (status == FastqPairReader::ReadStatus::END_OF_FILE) break;
            if (status == FastqPairReader::ReadStatus::READ_ERROR) {
                std::cerr << "Error: malformed or truncated FASTQ around pair\n";
                return 1;
            }

            ParsedBarcode cell_barcode = parse_barcodes_from_r1(pair.r1, cell_barcode_set);

            ParsedAntibody antibody_barcode = parse_antibody_from_r2(pair.r2, antibody_barcode_set);

            if (!cell_barcode.valid) {
                // std::cout << "\nPair #" << printed_pairs + 1 << " → no valid barcode\n";
            } else {
                // std::cout << "\nPair #" << printed_pairs + 1 << '\n';
                // std::cout << "Barcodes: bc1=" << cell_barcode.bc1 << " bc2=" << cell_barcode.bc2 << "\n";
                num_with_barcodes++;
            }

            if (!antibody_barcode.valid) {
                // std::cout << "Antibody barcode: NONE / INVALID \n";
            } else {
                // std::cout << "Antibody barcode (canonical 15 bp)\n" << antibody_barcode.barcode << "\n";
                num_with_ab_payload++;
            }

            if (cell_barcode.valid && antibody_barcode.valid) {
                std::string cell_id = cell_barcode.bc1 + "_" + cell_barcode.bc2;
                // std::cout << "Cell ID: " << cell_id << "\n";
                num_with_both++;
                counts[cell_id][antibody_barcode.barcode] += 1;
            }

            // std::cout << pair << "\n";

            printed_pairs++;
        }
        std::cout << "\n Summary over " << printed_pairs << " pairs:\n";
        std::cout << num_with_barcodes << " have valid cell barcodes\n"; // shows measurable improvement with hamming dictionary
        std::cout << num_with_ab_payload << " have valid ab payloads\n";
        std::cout << num_with_both << " have both cell + ab\n";

        
        const std::size_t MIN_COUNT = 10;
        std::string out_dir = "cell_counts";
        fs::create_directories(out_dir);

        // Print through counts in map.
        std::cout << "\n cell_id     antibody_barcode     count\n";
        for (const auto &[cell_id, antibody_counts] : counts) {

            std::size_t kept_antibodies = 0;
            for (const auto &[antibody_barcode, count] : antibody_counts) {
                if (count >= MIN_COUNT) kept_antibodies++;

                // std::cout << cell_id << ','
                //           << antibody_name << ','
                //           << count << '\n';
            }

            // If no antibodies kept, skip making a file.
            if (kept_antibodies == 0) continue;

            std::string filename = "cell_" + cell_id + ".csv";
            fs::path filepath = fs::path(out_dir) / filename;
            std::ofstream out(filepath);

            if (!out) {
                std::cerr << "Warning: could not open " << filepath << "for writing \n";
                continue;
            }

            out << "antibody_name, count\n";

            for (const auto &[antibody_barcode, count] : antibody_counts) {
                if (count < MIN_COUNT) continue;
                const std::string &antibody_name = antibody_barcode_to_name.at(antibody_barcode);
                out << antibody_name << "," << count << "\n";
            }
        }
        
    }
    catch (const std::exception &e)
    {
        std::cerr << "Fatal error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}