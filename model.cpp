#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <limits>
#include <random>
#include <stdexcept>
#include <cstdint>

std::vector<int> ComputeMinHashSignature(const std::vector<int>& expressed_genes, int num_hashes) {
    if (num_hashes <= 0) {
        throw std::invalid_argument("num_hashes must be > 0");
    }

    const int prime = 2147483647; // prime (2^31 - 1)
    std::vector<int> signature(num_hashes, std::numeric_limits<int>::max());

    // We want to reuse the coefficents to ensure the same hash functions are applied across different bins
    static std::vector<int> a_coeffs;
    static std::vector<int> b_coeffs;

    if (static_cast<int>(a_coeffs.size()) < num_hashes) {
        std::mt19937 rng(42); // fixed seed for reproducibility
        std::uniform_int_distribution<int> dist_a(1, prime - 1);
        std::uniform_int_distribution<int> dist_b(0, prime - 1);

        a_coeffs.resize(num_hashes);
        b_coeffs.resize(num_hashes);

        for (int i = 0; i < num_hashes; ++i) {
            a_coeffs[i] = dist_a(rng);
            b_coeffs[i] = dist_b(rng);
        }
    }

    // Standard Approach: Generating random coefficients A and B for using
    // h(x) = (A * x + B) % Prime_Modulo

    for (int gene_id : expressed_genes) {
        int x = ((gene_id % prime) + prime) % prime; // normalize gene_id

        for (int i = 0; i < num_hashes; ++i) {
            std::int64_t hashed = (static_cast<std::int64_t>(a_coeffs[i]) * x + b_coeffs[i]) % prime;
            if (hashed < signature[i]) {
                signature[i] = static_cast<int>(hashed);
            }
        }
    }
    return signature;
}

/** Function take in 2 inputs: 2 MinHash signatures
 *  It outputs the Jaccard similarity between the two signatures, i.e. how similar are the sigs
 *  Output in [0,1], 0 means identical, 1 means disjoint. 
 */
double SignatureSimilarity(const std::vector<int>& sig1, const std::vector<int>& sig2) {
    if (sig1.size() != sig2.size() || sig1.empty()) return 0.0;
    int matches = 0;
    for (size_t i = 0; i < sig1.size(); ++i) {
        if (sig1[i] == sig2[i]) ++matches;
    }
    return static_cast<double>(matches) / static_cast<double>(sig1.size());
}

/** Abstract base class for spatial sketching methods.
 *  We can add more functions here when we extend to othet sketching methods (spatially-aware variants for example)
 */
class SpatialSketcher {
public:
    virtual std::vector<std::size_t> ComputeSketch(int k) = 0;
    virtual ~SpatialSketcher() {} 
};

/** Baseline sketching implementation using MinHash signatures.
 *  Input: sparse_gene_matrix - N x G matrix (N bins, G number of expressed genes), spatial_coordinates_matrix - N x 2 matrix (N bins, 2D spatial coordinates
 *  Output: vector of indices of the selected bins in the sketch)
 * 
 *  Idea: First we compute the MinHash signatures for all the bins. These are just compact representations of the gene expression profiles. 
 *  Starting with a baseline bin (bin0), we select the bin that is farther from it (most dissimilar). This is measured according to the Jaccard similarity. 
 *  After we store the selected bin. Importantly, we also store the minimum distance of all unselected bins to any of the selected bins. 
 *  This absultely essential since we avoid N*k^2 comparisons by not computing distance to all selected bins everytime (only for the bin last added).
 */
class BaselineSketching : public SpatialSketcher {
private:
    std::vector<std::vector<int>> sparse_gene_matrix; // only stores indices of non-zero entries (expressed genees) for each bin
    std::vector<std::vector<double>> spatial_coordinates_matrix; 

public: 
    // COnstructor
    BaselineSketching(const std::vector<std::vector<int>>& gene_matrix, const std::vector<std::vector<double>>& spatial_coords)
        : sparse_gene_matrix(gene_matrix), spatial_coordinates_matrix(spatial_coords) {}

    /** Selects k bins that are most dissimilar in terms of gene expression profiles, using MinHash signatures and farthest-first selection. */
    std::vector<std::size_t> ComputeSketch(int k) override {
        std::cout << "Computing baseline MinHash sketch..." << std::endl;
        if (k <= 0 || sparse_gene_matrix.empty()) {
            return {};
        }

        const std::size_t total_bins = sparse_gene_matrix.size();
        const std::size_t target_k = static_cast<std::size_t>(k) < total_bins ? static_cast<std::size_t>(k) : total_bins;
        
        int num_hashes = 128; 
        std::vector<std::vector<int>> all_sigs;
        all_sigs.reserve(total_bins);

        // Compute our MinHash signatures for all bins upfront
        for (const auto& bin_genes : sparse_gene_matrix) {
            all_sigs.push_back(ComputeMinHashSignature(bin_genes, num_hashes));
        }

        //________________
        // Farthest-first selection on MinHash distance: distance = 1 - similarity.
        std::vector<std::size_t> sketch_indices;
        sketch_indices.reserve(target_k);

        std::vector<bool> selected(total_bins, false);
        std::vector<double> min_distance_to_selected(total_bins, 0.0);  // This will store the minimum distance of each unselected bin to any selected bin. We update this iteratively as we add new bins to the sketch.

        // Baseline: 
        const std::size_t seed_ind = 0;
        sketch_indices.push_back(seed_ind);
        selected[seed_ind] = true;

        for (std::size_t i = 0; i < total_bins; ++i) {
            if (selected[i]) continue;
            min_distance_to_selected[i] = 1.0 - SignatureSimilarity(all_sigs[i], all_sigs[seed_ind]);
        }

        // General Case: Select the bin that is farthest from any selected bin
        while (sketch_indices.size() < target_k) {
            std::size_t best_index = total_bins;
            double best_distance = -1.0;

            // Scan all unselected bins and pick the one whose nearest selected neighbor
            for (std::size_t i = 0; i < total_bins; ++i) {
                if (selected[i]) continue;
                if (min_distance_to_selected[i] > best_distance) {
                    best_distance = min_distance_to_selected[i];
                    best_index = i;
                }
            }
            // Stoping Condition
            if (best_index == total_bins) break;
            sketch_indices.push_back(best_index);
            selected[best_index] = true;

            // After adding best_index to the sketch, update each remaining bin's nearest-selected distance. We only compare to the newly added bin, then keep the minimum distance seen so far.
            for (std::size_t i = 0; i < total_bins; ++i) {
                if (selected[i]) continue;
                const double distance_to_new = 1.0 - SignatureSimilarity(all_sigs[i], all_sigs[best_index]);
                if (distance_to_new < min_distance_to_selected[i]) {
                    min_distance_to_selected[i] = distance_to_new;
                }
            }
        }

        return sketch_indices;
    }
};