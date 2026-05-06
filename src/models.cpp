#include "models.hpp"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <limits>
#include <random>
#include <stdexcept>
#include <cstdint>

// ________________________________________________________________
// _______________________Helper Functions_______________________ 

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

// ________________________________________________________________
// _______________________Baseline Sketching_______________________ 

BaselineSketching::BaselineSketching(const std::vector<std::vector<int>>& gene_matrix, const std::vector<std::vector<double>>& spatial_coords)
    : sparse_gene_matrix(gene_matrix), spatial_coordinates_matrix(spatial_coords) {}

std::vector<std::size_t> BaselineSketching::ComputeSketch(int k) {
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

    // Farthest-first selection on MinHash distance
    std::vector<std::size_t> sketch_indices;
    sketch_indices.reserve(target_k);

    std::vector<bool> selected(total_bins, false);
    std::vector<double> min_distance_to_selected(total_bins, 0.0);  

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
        
        // Stopping Condition
        if (best_index == total_bins) break;
        sketch_indices.push_back(best_index);
        selected[best_index] = true;

        // Update the minimum distance for remaining bins relative to the newly added bin
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

// ________________________________________________________________
// _________________________________QuadTree Implementation_________________________________

bool BoundingBox::Contains(double x, double y) const {
    return (x >= x_min && x <= x_max && y >= y_min && y <= y_max);
}

QuadTreeNode::QuadTreeNode(BoundingBox b, int cap) : boundary(b), capacity(cap), 
    northWest(nullptr), northEast(nullptr), southWest(nullptr), southEast(nullptr) {}

QuadTreeNode::~QuadTreeNode() {
    delete northWest;
    delete northEast;
    delete southWest;
    delete southEast;
}

void QuadTreeNode::Subdivide() {
    double x_mid = boundary.x_min + (boundary.x_max - boundary.x_min) / 2.0;
    double y_mid = boundary.y_min + (boundary.y_max - boundary.y_min) / 2.0;

    northWest = new QuadTreeNode({boundary.x_min, x_mid, boundary.y_min, y_mid}, capacity);
    northEast = new QuadTreeNode({x_mid, boundary.x_max, boundary.y_min, y_mid}, capacity);
    southWest = new QuadTreeNode({boundary.x_min, x_mid, y_mid, boundary.y_max}, capacity);
    southEast = new QuadTreeNode({x_mid, boundary.x_max, y_mid, boundary.y_max}, capacity);
}

bool QuadTreeNode::Insert(std::size_t index, double x, double y) {
    if (!boundary.Contains(x, y)) {
        return false;
    }

    if (northWest == nullptr && points.size() < static_cast<std::size_t>(capacity)) {
        points.push_back({index, x, y}); 
        return true;
    }

    // If the node is full, subdivide it
    // FIX: If we subdivide, we MUST push the points from the parent node down into the children
    if (northWest == nullptr) {
        Subdivide();
        
        for (const auto& p : points) {
            if (northWest->Insert(p.index, p.x, p.y)) continue;
            if (northEast->Insert(p.index, p.x, p.y)) continue;
            if (southWest->Insert(p.index, p.x, p.y)) continue;
            southEast->Insert(p.index, p.x, p.y);
        }
        // Clear the parent node's points so data only lives in leaves
        points.clear(); 
    }
    
    if (northWest->Insert(index, x, y)) return true;
    if (northEast->Insert(index, x, y)) return true;
    if (southWest->Insert(index, x, y)) return true;
    if (southEast->Insert(index, x, y)) return true;

    return false;
}

void QuadTreeSketching::AssignLeafBudgets(QuadTreeNode* node, std::vector<int>& bin_to_leaf, std::vector<int>& leaf_budgets, int max_per_leaf) {
    if (!node) return;

    // Check for nullptr (leaf node)
    if (!node->northWest) {
        int leaf_id = leaf_budgets.size();
        leaf_budgets.push_back(max_per_leaf);
        
        // Record which leaf each bin belongs to, so we can decrease the budget as we go
        for (const auto& p : node->points) { 
            bin_to_leaf[p.index] = leaf_id;  
        }
        return;
    }
    
    AssignLeafBudgets(node->northWest, bin_to_leaf, leaf_budgets, max_per_leaf);
    AssignLeafBudgets(node->northEast, bin_to_leaf, leaf_budgets, max_per_leaf);
    AssignLeafBudgets(node->southWest, bin_to_leaf, leaf_budgets, max_per_leaf);
    AssignLeafBudgets(node->southEast, bin_to_leaf, leaf_budgets, max_per_leaf);
}

QuadTreeSketching::QuadTreeSketching(const std::vector<std::vector<int>>& gene_matrix, 
    const std::vector<std::vector<double>>& spatial_coords, int capacity, int budget)
    : sparse_gene_matrix(gene_matrix), spatial_coordinates_matrix(spatial_coords), node_capacity(capacity), leaf_budget(budget) {}

std::vector<std::size_t> QuadTreeSketching::ComputeSketch(int k) {
    std::cout << "Computing Spatially-Constrained QuadTree Sketch..." << std::endl;
    if (k <= 0 || sparse_gene_matrix.empty()) return {};

    const std::size_t total_bins = sparse_gene_matrix.size();
    const std::size_t target_k = static_cast<std::size_t>(k) < total_bins ? static_cast<std::size_t>(k) : total_bins;
    
    std::cout << "  -> Building physical spatial tree..." << std::endl;
    
    // Find the boundaries of the tissue slide
    double min_x = std::numeric_limits<double>::max();
    double max_x = std::numeric_limits<double>::lowest();
    double min_y = std::numeric_limits<double>::max();
    double max_y = std::numeric_limits<double>::lowest();

    for (const auto& coords : spatial_coordinates_matrix) {
        if (coords[0] < min_x) min_x = coords[0];
        if (coords[0] > max_x) max_x = coords[0];
        if (coords[1] < min_y) min_y = coords[1];
        if (coords[1] > max_y) max_y = coords[1];
    }

    // UPDATE: Added capacity as paramter
    int final_capacity;
    if (node_capacity > 0) {
        final_capacity = node_capacity;
        std::cout << "     User Defined Capacity: " << final_capacity << std::endl;
    } else {
        final_capacity = (total_bins / target_k) + 1; 
        std::cout << "     Dynamic Capacity: " << final_capacity << std::endl;
    }

    QuadTreeNode* root = new QuadTreeNode({min_x, max_x, min_y, max_y}, final_capacity);

    // Add all 480,000 bins into the tree 
    for (std::size_t i = 0; i < total_bins; ++i) {
        root->Insert(i, spatial_coordinates_matrix[i][0], spatial_coordinates_matrix[i][1]);
    }

    std::vector<int> bin_to_leaf_id(total_bins, -1);
    std::vector<int> leaf_budgets;

    int final_budget;
    
    // UPDATED: Added budget as a parameter
    if (leaf_budget > 0) {
        final_budget = leaf_budget;
        std::cout << "     User Defined Budget: " << final_budget << std::endl;
    } else {
        int estimated_leaves = total_bins / final_capacity;
        if (estimated_leaves == 0) estimated_leaves = 1;
        final_budget = (target_k / estimated_leaves) + 2; 
        std::cout << "     Dynamic Budget: " << final_budget << std::endl;
    }
    
    // Cap each leaf at 2 selectable points 
    AssignLeafBudgets(root, bin_to_leaf_id, leaf_budgets, final_budget); 

        //_______________________
    // THE MINHASH SEARCH
    std::cout << "  -> Running transcriptomic search..." << std::endl;
    
    int num_hashes = 128; 
    std::vector<std::vector<int>> all_sigs;
    all_sigs.reserve(total_bins);

    for (const auto& bin_genes : sparse_gene_matrix) {
        all_sigs.push_back(ComputeMinHashSignature(bin_genes, num_hashes));
    }

    std::vector<std::size_t> sketch_indices;
    sketch_indices.reserve(target_k);
    std::vector<bool> selected(total_bins, false);
    std::vector<double> min_distance_to_selected(total_bins, 0.0);

    // Seed initialization
    const std::size_t seed_ind = 0;
    sketch_indices.push_back(seed_ind);
    selected[seed_ind] = true;
    leaf_budgets[bin_to_leaf_id[seed_ind]]--; // Deduct from the seed's spatial leaf budget!

    for (std::size_t i = 0; i < total_bins; ++i) {
        if (selected[i]) continue;
        min_distance_to_selected[i] = 1.0 - SignatureSimilarity(all_sigs[i], all_sigs[seed_ind]);
    }

    // General Case: Farthest-first search
    while (sketch_indices.size() < target_k) {
        std::size_t best_index = total_bins;
        double best_distance = -1.0;

        for (std::size_t i = 0; i < total_bins; ++i) {
            if (selected[i]) continue;
            
            // Ignore block if oversampled
            if (leaf_budgets[bin_to_leaf_id[i]] <= 0) continue;

            if (min_distance_to_selected[i] > best_distance) {
                best_distance = min_distance_to_selected[i];
                best_index = i;
            }
        }

        if (best_index == total_bins) {
            std::cerr << "Warning: Ran out of spatial budget before reaching k points." << std::endl;
            break;
        }

        sketch_indices.push_back(best_index);
        selected[best_index] = true;
        leaf_budgets[bin_to_leaf_id[best_index]]--; // Deduct the point from this leaf's budget

        // Update distances for remaining bins
        for (std::size_t i = 0; i < total_bins; ++i) {
            if (selected[i]) continue;
            const double distance_to_new = 1.0 - SignatureSimilarity(all_sigs[i], all_sigs[best_index]);
            if (distance_to_new < min_distance_to_selected[i]) {
                min_distance_to_selected[i] = distance_to_new;
            }
        }
    }

    delete root;
    return sketch_indices;
}