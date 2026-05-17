#pragma once

#include <vector>
#include <cstddef>

// ________________________________________________________________
// _______________________Evaluation Metrics_______________________

/** Helper function to extract a specific quantile from a vector of values */
double GetQuantile(std::vector<double>& values, double quantile);

/** Computes the Robust Transcriptomic Hausdorff Distance (95th percentile)
 * using MinHash Jaccard approximation. Lower is better. */
double ComputeTranscriptomicHausdorff(
    const std::vector<std::vector<int>>& sparse_gene_matrix,
    const std::vector<std::size_t>& sketch_indices,
    int num_hashes = 128,
    double quantile = 0.95);

/** Computes the Robust Coordinate Hausdorff Distance (95th percentile)
 * using Euclidean physical distance. Lower is better. */
double ComputeCoordinateHausdorff(
    const std::vector<std::vector<double>>& spatial_coords,
    const std::vector<std::size_t>& sketch_indices,
    double quantile = 0.95);

// ________________________________________________________________
// _______________________MinHash Helpers__________________________________

std::vector<int> ComputeMinHashSignature(const std::vector<int>& expressed_genes, int num_hashes);

/** Function take in 2 inputs: 2 MinHash signatures
 *  It outputs the Jaccard similarity between the two signatures, i.e. how similar are the sigs
 *  Output in [0,1], 0 means identical, 1 means disjoint. 
 */
double SignatureSimilarity(const std::vector<int>& sig1, const std::vector<int>& sig2);

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
    std::vector<std::vector<int>> sparse_gene_matrix; 
    std::vector<std::vector<double>> spatial_coordinates_matrix; 

public: 
    BaselineSketching(const std::vector<std::vector<int>>& gene_matrix, const std::vector<std::vector<double>>& spatial_coords);
    
    /** Selects k bins that are most dissimilar in terms of gene expression profiles, using MinHash signatures and farthest-first selection. */
    std::vector<std::size_t> ComputeSketch(int k) override;
};

// _________________________________QuadTree Implementation_________________________________
// Idea: Want to enforce uniform sampling across the spatial domain

// Defines the 2D physical boundaries of any given node
struct BoundingBox {
    double x_min, x_max;
    double y_min, y_max;

    bool Contains(double x, double y) const;
};

struct SpatialPoint {
    std::size_t index;
    double x;
    double y;
};

class QuadTreeNode {
public:
    BoundingBox boundary;
    std::vector<SpatialPoint> points;   // stores (index, x, y) for points in this node
    int capacity; 

    QuadTreeNode* northWest;
    QuadTreeNode* northEast;
    QuadTreeNode* southWest;
    QuadTreeNode* southEast;

    QuadTreeNode(BoundingBox b, int cap);
    ~QuadTreeNode();

    void Subdivide();
    bool Insert(std::size_t index, double x, double y);
};

// QuadTree sketching implementation
class QuadTreeSketching : public SpatialSketcher {
private:
    std::vector<std::vector<int>> sparse_gene_matrix; 
    std::vector<std::vector<double>> spatial_coordinates_matrix; 
    int node_capacity;
    int leaf_budget;

    void AssignLeafBudgets(QuadTreeNode* node, std::vector<int>& bin_to_leaf, std::vector<int>& leaf_budgets, int max_per_leaf);

public: 
    QuadTreeSketching(const std::vector<std::vector<int>>& gene_matrix, const std::vector<std::vector<double>>& spatial_coords,
    int capacity=-1, int leaf_budget=-1);
    std::vector<std::size_t> ComputeSketch(int k) override;
};

// Add to models.hpp

struct VoronoiCell {
    double seed_x, seed_y;
    int cell_id;
};

class VoronoiSketching : public SpatialSketcher {
public:
    // Default args belong ONLY in the header, never in the .cpp definition
    VoronoiSketching(const std::vector<std::vector<int>>& gene_matrix,
                     const std::vector<std::vector<double>>& spatial_coords,
                     int num_cells = -1,
                     double epsilon = 0.1);

    std::vector<std::size_t> ComputeSketch(int k) override;

private:
    const std::vector<std::vector<int>>& sparse_gene_matrix;
    const std::vector<std::vector<double>>& spatial_coordinates_matrix;
    int num_voronoi_cells;
    double epsilon;

    std::vector<VoronoiCell> GenerateSeeds(double min_x, double max_x,
                                           double min_y, double max_y);
    int AssignToCell(double x, double y, const std::vector<VoronoiCell>& seeds);
};