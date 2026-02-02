#ifndef BARCODE_INDEX_H
#define BARCODE_INDEX_H

#include <string>
#include <unordered_set>
#include <unordered_map>

class BarcodeIndex {

public:
    // Prevents implicit conversion.
    explicit BarcodeIndex(const std::string& csv_path);

    bool is_valid(const std::string& bc) const; // const at the end here means this method is read-only. Doesn't impact class members.

    bool find_canonical_barcode(const std::string& observed, std::string& canonical) const; // read only method, doesn't impact class members.

    // Get number of canonical barcodes loaded
    std::size_t size() const { return _canonical_barcodes_.size(); }
    std::size_t hamming_dict_size() const { return _noisy_to_canonical_.size(); }

private:
    std::unordered_map <std::string, std::string> _noisy_to_canonical_;
    std::unordered_set <std::string> _canonical_barcodes_;
    void add_hamming_neighbors(const std::string& bc, int hamming_dist);
};

#endif // BARCODE_INDEX_H