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
        const int max_pairs = 1000000;
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

        std::cout << "\nSummary over " << printed_pairs << " pairs:\n";
        std::cout << num_with_barcodes << " have valid cell barcodes\n";
        std::cout << num_with_ab_payload << " have valid ab payloads\n";
        std::cout << num_with_both << " have both cell + ab\n";

        // Output single CSV file
        std::string output_file = "antibody_counts_orig.csv";
        std::ofstream out(output_file);

        if (!out) {
            std::cerr << "Error: could not open " << output_file << " for writing\n";
            return 1;
        }

        // Write header
        out << "cell_id,antibody_name,count\n";

        std::size_t total_rows = 0;

        const std::size_t MIN_COUNT = 10;

        for (const auto &[cell_id, antibody_counts] : counts) {
            for (const auto &[antibody_barcode, count] : antibody_counts) {
                if (count < MIN_COUNT) continue;
                
                auto it = antibody_barcode_to_name.find(antibody_barcode);
                std::string antibody_name = (it != antibody_barcode_to_name.end()) ? it->second : "UNKNOWN";

                out << cell_id << "," << antibody_name << "," << count << "\n";
                total_rows++;
            }
        }

        out.close();

        std::cout << "\nWrote " << total_rows << " rows to " << output_file << "\n";
        std::cout << "Unique cells: " << counts.size() << "\n";
    }
    catch (const std::exception &e)
    {
        std::cerr << "Fatal error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}