#include <iostream>
#include "fastq_reader.h"
#include "barcode_index.h"
#include "dabseq_utilities.h"
#include <fstream>
#include <iomanip>    // for std::setprecision, std::setw, std::left
#include <algorithm>  // for std::sort, std::min
#include <vector>     // for std::vector

int main(int argc, char **argv)
{
    std::cout << "========================================\n";
    std::cout << "  DAb-seq C++ Pipeline\n";
    std::cout << "========================================\n\n";

    if (argc != 5)
    {
        std::cerr << "Usage: " << argv[0] << " R1.fastq[.gz] R2.fastq[.gz] cell_barcodes.csv antibody_barcodes.csv\n";
        return 1;
    }

    const std::string r1_path = argv[1];
    const std::string r2_path = argv[2];
    const std::string cell_barcodes_csv = argv[3];
    const std::string antibody_barcodes_csv = argv[4];


    std::cout << "[Input Files]\n";
    std::cout << "  R1 FASTQ:           " << r1_path << "\n";
    std::cout << "  R2 FASTQ:           " << r2_path << "\n";
    std::cout << "  Cell barcodes:      " << cell_barcodes_csv << "\n";
    std::cout << "  Antibody barcodes:  " << antibody_barcodes_csv << "\n\n";

    try
    {
        // Load cell barcode whitelist.
        std::cout << "[Loading Barcodes]\n";
        std::cout << "  Loading cell barcodes..." << std::flush;
        BarcodeIndex cell_barcode_set(cell_barcodes_csv);
        std::cout << " done.\n";
        std::cout << "    -> " << cell_barcode_set.size() << " canonical barcodes\n";
        std::cout << "    -> " << cell_barcode_set.hamming_dict_size() << " entries in hamming dictionary\n";


        // Load antibody barcode whitelist.
        std::cout << "  Loading antibody barcodes..." << std::flush;
        BarcodeIndex antibody_barcode_set(antibody_barcodes_csv);
        std::cout << " done.\n";
        std::cout << "    -> " << antibody_barcode_set.size() << " canonical barcodes\n";
        std::cout << "    -> " << antibody_barcode_set.hamming_dict_size() << " entries in hamming dictionary\n";

        std::cout << "  Loading antibody name map..." << std::flush;
        std::unordered_map<std::string, std::string> antibody_barcode_to_name = load_antibody_name_map(antibody_barcodes_csv);
        std::cout << " done.\n";
        std::cout << "    -> " << antibody_barcode_to_name.size() << " antibody name mappings\n\n";

        std::cout << "[Opening FASTQ Files]\n";
        FastqPairReader reader(r1_path, r2_path);
        std::cout << "  FASTQ files opened successfully.\n\n";

        FastqPairReader::FastqPair pair;

        std::size_t total_pairs = 0;
        const int max_pairs = 1000000;
        std::size_t num_with_barcodes = 0; // measure effectiveness of hamming dictionary
        std::size_t num_with_ab_payload = 0;
        std::size_t num_with_both = 0;

        //                  cell_id     antibody, cnt
        std::unordered_map<std::string, std::unordered_map<std::string, std::size_t>> counts;


        std::cout << "[Processing Reads]\n";
        std::cout << "  Processing";
        if (max_pairs > 0) {
            std::cout << " (max " << max_pairs << " pairs)";
        }
        std::cout << "...\n";

        // Progress tracking
        const std::size_t progress_interval = 1000000;

        while (max_pairs == 0 || total_pairs < max_pairs) // check logic
        {
            FastqPairReader::ReadStatus status = reader.next_record(pair);

            if (status == FastqPairReader::ReadStatus::END_OF_FILE) {
                std::cout << "  Reached end of file.\n";
                break;
            }
            
            if (status == FastqPairReader::ReadStatus::READ_ERROR) {
                std::cerr << "\n  Error: malformed or truncated FASTQ at pair " << total_pairs + 1 << "\n"; // because haven't incremented total pairs yet
                return 1;
            }

            total_pairs++;

            // Progress indicator
            if (total_pairs % progress_interval == 0) {
                std::cout << "  Processed " << total_pairs << " pairs...\r" << std::flush;
            }

            // Parse cell barcode from R1
            ParsedBarcode cell_barcode = parse_barcodes_from_r1(pair.r1, cell_barcode_set);
            
            // Parse antibody barcode from R2
            ParsedAntibody antibody_barcode = parse_antibody_from_r2(pair.r2, antibody_barcode_set);

            if (cell_barcode.valid) {
                num_with_barcodes++;
            }

            if (antibody_barcode.valid) {
                num_with_ab_payload++;
            }

            if (cell_barcode.valid && antibody_barcode.valid) {

                std::string cell_id = cell_barcode.bc1 + "_" + cell_barcode.bc2;
                num_with_both++;
                counts[cell_id][antibody_barcode.barcode] += 1;
            }
        }

        // Clear progress line and print final count
        std::cout << "  Processed " << total_pairs << " pairs total.        \n\n";

        // Summary statistics
        std::cout << "[Summary Statistics]\n";
        std::cout << "  Total read pairs processed:    " << total_pairs << "\n";
        std::cout << "  Valid cell barcodes:           " << num_with_barcodes 
                  << " (" << std::fixed << std::setprecision(1) 
                  << (100.0 * num_with_barcodes / total_pairs) << "%)\n";
        std::cout << "  Valid antibody payloads:       " << num_with_ab_payload
                  << " (" << std::fixed << std::setprecision(1)
                  << (100.0 * num_with_ab_payload / total_pairs) << "%)\n";
        std::cout << "  Both valid (countable reads):  " << num_with_both
                  << " (" << std::fixed << std::setprecision(1)
                  << (100.0 * num_with_both / total_pairs) << "%)\n";
        std::cout << "  Unique cell barcodes observed: " << counts.size() << "\n\n";


// Output file
        std::string output_file = "antibody_counts.tsv";

        std::cout << "[Writing Output File]\n";
        std::cout << "  Output file: " << output_file << "\n";

        std::ofstream out(output_file);
        if (!out) {
            std::cerr << "  Error: could not open " << output_file << " for writing\n";
            return 1;
        }

        // Write header
        out << "cell_id\tcell_bc1\tcell_bc2\tantibody_barcode\tantibody_name\tcount\n";

        std::size_t total_rows = 0;
        std::size_t cells_written = 0;

        // Sort cells by ID for consistent output
        std::vector<std::string> sorted_cell_ids;
        for (const auto &[cell_id, _] : counts) {
            sorted_cell_ids.push_back(cell_id);
        }
        std::sort(sorted_cell_ids.begin(), sorted_cell_ids.end());

        for (const auto &cell_id : sorted_cell_ids) {
            const auto &antibody_counts = counts.at(cell_id);
            
            // Parse bc1 and bc2 from cell_id (format: "BC1_BC2")
            std::size_t underscore_pos = cell_id.find('_');
            std::string bc1 = (underscore_pos != std::string::npos) ? cell_id.substr(0, underscore_pos) : cell_id;
            std::string bc2 = (underscore_pos != std::string::npos) ? cell_id.substr(underscore_pos + 1) : "";

            bool cell_has_output = false;

            // Sort antibodies by count descending for this cell
            std::vector<std::pair<std::string, std::size_t>> sorted_abs(antibody_counts.begin(), antibody_counts.end());
            std::sort(sorted_abs.begin(), sorted_abs.end(),
                      [](const auto& a, const auto& b) { return a.second > b.second; });

            for (const auto &[antibody_barcode, count] : sorted_abs) {
                if (count < 10) continue;
                // Look up antibody name
                std::string antibody_name;
                auto it = antibody_barcode_to_name.find(antibody_barcode);
                if (it != antibody_barcode_to_name.end()) {
                    antibody_name = it->second;
                } else {
                    antibody_name = "UNKNOWN";
                }

                out << cell_id << "\t" 
                    << bc1 << "\t" 
                    << bc2 << "\t"
                    << antibody_barcode << "\t" 
                    << antibody_name << "\t" 
                    << count << "\n";
                
                total_rows++;
                cell_has_output = true;
            }

            if (cell_has_output) {
                cells_written++;
            }
        }

        out.close();

        std::cout << "  Cells written:     " << cells_written << "\n";
        std::cout << "  Total rows:        " << total_rows << "\n\n";

        // Print top cells by total counts (optional detailed output)
        std::cout << "[Top Cells by Total Counts]\n";
        
        // Create vector for sorting
        std::vector<std::pair<std::string, std::size_t>> cell_totals;
        for (const auto &[cell_id, antibody_counts] : counts) {
            std::size_t total = 0;
            for (const auto &[ab, cnt] : antibody_counts) {
                total += cnt;
            }
            cell_totals.emplace_back(cell_id, total);
        }
        
        // Sort by count descending
        std::sort(cell_totals.begin(), cell_totals.end(),
                  [](const auto& a, const auto& b) { return a.second > b.second; });
        
        // Print top 10
        std::cout << "  " << std::left << std::setw(25) << "Cell ID" << "Total Counts\n";
        std::cout << "  " << std::string(40, '-') << "\n";
        for (std::size_t i = 0; i < std::min<std::size_t>(10, cell_totals.size()); i++) {
            std::cout << "  " << std::left << std::setw(25) << cell_totals[i].first 
                      << cell_totals[i].second << "\n";
        }

        std::cout << "\n========================================\n";
        std::cout << "  Pipeline Complete!\n";
        std::cout << "========================================\n";
    }
    catch (const std::exception &e)
    {
        std::cerr << "\nFatal error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
