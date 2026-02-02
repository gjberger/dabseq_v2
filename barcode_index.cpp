#include "barcode_index.h"
#include <fstream>
#include <sstream>
#include <stdexcept>

/* Building a hamming dictionary to allow for hamming distance 1/2. Mission Bio
 * docs guarantee that the cell barcodes are more than 3 Levenshtein distance apart.
 * https://missionbio.com/wp-content/uploads/2019/10/WhitePaper_MissionBio_TapestriPlatform_RevA.pdf
 * So, we are safe in assuming there will be no collisions if we build a hash table of noisy barcodes
 * back to canonical barcodes.
 */


/**
 * @brief Constructor, Barcode Index object.
 * @param csv_path path to csv of antibody/cell barcodes.
 * Loads either:
 * -> 1536 Mission Bio Cell Barcodes (halfs)
 * -> 46 TotalSeq-B antibody barcodes
 */
BarcodeIndex::BarcodeIndex(const std::string &csv_path) {
    std::ifstream barcode_csv(csv_path); // Open CSV file.

    if (!barcode_csv) { // Throw error if CSV doesn't open.
        throw std::runtime_error("Failed to open barcode CSV: " + csv_path);
    }

    std::string line; // Temp. holder for each line.

    // Loop and parse lines of CSV.
    while (std::getline(barcode_csv, line)) {
        if (line.empty()) continue; // Skip blank lines.
        std::size_t comma_pos = line.find(','); // CSV comma-delimited, find comma.
        if (comma_pos == std::string::npos) { // Throw error if no comma.
            throw std::runtime_error("Malformed barcode line (no comma): " + line);
        }
        
        // CSV structure is:
        // -> (1) cell_bc, #
        // -> (2) antibody_bc, antibody_name
        std::string bc = line.substr(0, comma_pos); // Parse barcode from line.

        // ID in CSV list -> potentially build an unordered map.
        // std::string id_str = line.substr(comma_pos + 1);
        // int id = std::stoi(id_str); // assuming the ID is an integer, which it is.

        _canonical_barcodes_.insert(bc); // Insert barcode into unordered set.
        add_hamming_neighbors(bc, 1); // Building hamming neighbors, distance 1.
    }
}
/**
 * @brief check barcode validity.
 * Barcodes are stored in an unordered set, so O(nlogn) search time.
 * 
 * @param bc barcode.
 * @return true, barcode in whitelist
 * @return false, barcode not in whitelist.
 */
bool BarcodeIndex::is_valid(const std::string &bc) const {
    return _canonical_barcodes_.contains(bc);
}

/**
 * @brief check if barcode is in hamming neighbor list.
 * 
 * Cell/antibody barcodes are located at a constant area dependent on
 * Mission Bio chemistry. That area is extracted, and error-corrected
 * using the Hamming Neighbor List. This function maps an observed barcode
 * to it's canonical version. Barcodes are different enough that this causes
 * no collisions.
 * 
 * Eg:
 * -> Mission Bio Cell Barcode, antibody r1
 * Cell barcode location is immediately upstream of a read 1 start motif, "GTACTCGCAGTAGTC"
 * Take the following r1 read:
 * "TNGACCATGAGTACGTACGAGTCTGAACGGTTGTACTCGCAGTAGTCCGACTGAGATNCTAGATCGGCTCTAATGAGGAACACGGCCATGAGTNGCGTAGTCNTCTCGCNGNCTCTTATACACATCTCCGNGCCCACGANACAGGCAGAAA"
 * "TNGACCATGAGTACGTACGAGTCTGAACGGTT    GTACTCGCAGTAGTC   CGACTGAGATNCTAGATCGGCTCTAATGAGGAACACGGCCATGAGTNGCGTAGTCNTCTCGCNGNCTCTTATACACATCTCCGNGCCCACGANACAGGCAGAAA"
 *                                      ^^^^^^^^^^^^^^^
 * The motif begins at index 32 (0-indexed), meaning the barcode region is immediately upstream.
 * "TNGACCATGAGTACGTACGAGTCTGAACGGTT"
 * Part 1 of the barcode is the first 9 base pairs
 * Part 2 is the last 9 base pairs
 * "TNGACCATG   AGTACGTACGAGTC  TGAACGGTT"
 *  ^^^^^^^^^                   ^^^^^^^^^
 * Part 1 has N, which could be any base. The pipeline checks the hamming
 * dictionary and finds that cell barcode 10 matches, TAGACCATG.
 * 
 * Part 2 is a direct match with TGAACGGTT, barcode 506.
 * 
 * Thus, the full cell barcode is TAGACCATG_TGAACGGTT.
 * 
 * -> TotalSeq-b Antibody Barcode, antibody r2
 * The antibody barcodes are in read 2, and roughly structured like:
 * [junk] + [5' antibody handle] + [15 bp antibody barcode] + [more] + [3' antibody handle] // I think no more.
 * They are sandwiched immediately between the 5' and 3' handles.
 * The 5' handles can be either of "TGACTACGCTACTCATGG" or "TGACTACACTACTCATGG"
 * The 3' handle is "GAGCCGATCTAGTATCTCAGTCG"
 * Take the following r2 read:
 * "CGANATGACTACGCTACTCATGGCCGTGTTCCNCATTAGAGCCGATCTAGNATCTNAGTCNGACTACTNCGANTNCAANCGTNCAGNCTCGNANGNANTCNTNGTNTNCTNTNTCTTNTACACATNTGACNCTNCCNACGANANNGNTNGN"
 * "CGANA   TGACTACGCTACTCATGG  CCGTGTTCCNCATTAGAGCCGATCTAGNATCTNAGTCNGACTACTNCGANTNCAANCGTNCAGNCTCGNANGNANTCNTNGTNTNCTNTNTCTTNTACACATNTGACNCTNCCNACGANANNGNTNGN"
 *          ^^^^^^^^^^^^^^^^^^ <- 5' handle
 * 
 * The motif begins at index 5, so we focus on what comes after it:
 * "CCGTGTTCCNCATTA     GAGCCGATCTAG    NATCTNAGTCNGACTACTNCGANTNCAANCGTNCAGNCTCGNANGNANTCNTNGTNTNCTNTNTCTTNTACACATNTGACNCTNCCNACGANANNGNTNGN"
 *                      ^^^^^^^^^^^^ <- 3' handle prefix
 * The prefix of the 3' handle is a bit later, and the antibody handle is directly between them
 * CCGTGTTCCNCATTA
 * ^^^^^^^^^^^^^^^
 * This antibody barcode has an N which is resolved by the hamming dictionary
 * and maps to CCGTGTTCCTCATTA,CD71.
 * 
 * @param observed potentially noisy observed barcode
 * @param canonical corresponding canonical barcode, changed by reference
 * @return true if there is a matching hamming neighbor
 * @return false else
 */
bool BarcodeIndex::find_canonical_barcode(const std::string &observed, std::string &canonical) const {
    auto iterator = _noisy_to_canonical_.find(observed); // Find noisy barcode in map.
    if (iterator == _noisy_to_canonical_.end()) { // If it doesn't exist, ignore.
        return false;
    }
    canonical = iterator->second; // Canonical will be stored adjacent to observed.
    return true;
}

/**
 * @brief generates a map from noisy barcodes to canonical, if they exist.
 * 
 * @param bc barcode to generate hamming dictionary for.
 * @param hamming_dist allowed hamming dist.
 */
void BarcodeIndex::add_hamming_neighbors(const std::string &bc, int hamming_dist) {
    if (hamming_dist != 1) {
        throw std::runtime_error("Only hamming_dist = 1 supported right now");
    }
    static const char BASES[] = {'A', 'C', 'G', 'T', 'N'}; // static means there is only one instance of this in the whole program

    _noisy_to_canonical_.insert({bc, bc}); // Canonical barcode maps to itself.

    const std::size_t LENGTH = bc.size(); // Barcode length.

    for (std::size_t i = 0; i < LENGTH; i++) { // Loop over length of barcode, modifying each character.
        char original = bc[i];

        for (char base : BASES) { // Add a variant that is 1 base off for each base.
            if (base == original) continue; // skip over original, we only want different bases.

            std::string neighbor = bc; // Set neighbor to original barcode. --> is this 9 or 18?
            neighbor[i] = base; // Modify with different base.
            _noisy_to_canonical_.insert({neighbor, bc}); // Map noisy hamming neighbor to canonical barcode.
        }
    }
}