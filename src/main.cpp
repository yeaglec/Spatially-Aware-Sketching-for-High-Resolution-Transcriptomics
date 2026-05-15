#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <tuple>

#include "models.hpp"

// Helper to load the spatial coordinates CSV
std::vector<std::vector<double>> LoadSpatialCoords(const std::string& filename) {
    std::vector<std::vector<double>> coords;
    std::ifstream file(filename);
    std::string line, val;

    while (std::getline(file, line)) {
        std::vector<double> point;
        std::stringstream ss(line);
        while (std::getline(ss, val, ',')) {
            point.push_back(std::stod(val));
        }
        coords.push_back(point);
    }
    return coords;
}

// Helper to load the ragged sparse genes text file
std::vector<std::vector<int>> LoadSparseGenes(const std::string& filename) {
    std::vector<std::vector<int>> sparse_genes;
    std::ifstream file(filename);
    std::string line, val;

    while (std::getline(file, line)) {
        std::vector<int> bin_genes;
        if (!line.empty()) {
            std::stringstream ss(line);
            while (std::getline(ss, val, ',')) {
                bin_genes.push_back(std::stoi(val));
            }
        }
        sparse_genes.push_back(bin_genes);
    }
    return sparse_genes;
}

// CSV EXPORT HELPER
std::tuple<int, int, std::string> ExportToyMarkers(const std::vector<int>& genes) {
    int epithelial_markers = 0;
    int stromal_markers = 0;

    for (int g : genes) {
        if (g >= 1 && g <= 3) {
            ++epithelial_markers;
        } else if (g >= 4 && g <= 6) {
            ++stromal_markers;
        }
    }

    std::string dominant_program = "Mixed";
    if (epithelial_markers > stromal_markers) {
        dominant_program = "Epithelial";
    } else if (stromal_markers > epithelial_markers) {
        dominant_program = "Stromal";
    }

    return {epithelial_markers, stromal_markers, dominant_program};
}

// ______________________Toy Data Helper Functions______________________

void ExportToyResultsToCSV(const std::string& filename, 
                           const std::vector<std::vector<int>>& sparse_genes,
                           const std::vector<std::vector<double>>& spatial_coords,
                           const std::vector<std::size_t>& sketch_indices) {
    std::ofstream out_file(filename);
    if (!out_file.is_open()) return;

    std::vector<bool> is_sketched(spatial_coords.size(), false);
    for (std::size_t idx : sketch_indices) is_sketched[idx] = true;

    out_file << "Original_Index,X_Coordinate,Y_Coordinate,Epithelial_Markers,Stromal_Markers,Dominant_Program,Is_Sketched\n";

    for (std::size_t i = 0; i < spatial_coords.size(); ++i) {
        const auto [epi_count, stroma_count, dominant_program] = ExportToyMarkers(sparse_genes[i]);
        out_file << i << "," << spatial_coords[i][0] << "," << spatial_coords[i][1] << "," 
                 << epi_count << "," << stroma_count << "," << dominant_program << ","
                 << (is_sketched[i] ? 1 : 0) << "\n";
    }
    out_file.close();
}

void RunToyPipeline(const std::string& experiment_name, const std::vector<std::vector<int>>& sparse_genes, const std::vector<std::vector<double>>& spatial_coords, int k) {
    SpatialSketcher* my_sketcher = new BaselineSketching(sparse_genes, spatial_coords);
    std::vector<std::size_t> final_sketch = my_sketcher->ComputeSketch(k);
    ExportToyResultsToCSV(experiment_name + "_results.csv", sparse_genes, spatial_coords, final_sketch);
    delete my_sketcher;
}

// ______________________Visium Data Helper Functions______________________

void ExportVisiumResultsToCSV(const std::string& filename, 
                              const std::vector<std::vector<double>>& spatial_coords,
                              const std::vector<std::size_t>& sketch_indices) {
    
    std::ofstream out_file(filename);
    if (!out_file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << " for writing.\n";
        return;
    }

    // Create a boolean mask to tag which points were sketched
    std::vector<bool> is_sketched(spatial_coords.size(), false);
    for (std::size_t idx : sketch_indices) {
        is_sketched[idx] = true;
    }

    // Clean header matching the Visium visualization script
    out_file << "Original_Index,X_Coordinate,Y_Coordinate,Is_Sketched\n";

    // Write all data
    for (std::size_t i = 0; i < spatial_coords.size(); ++i) {
        out_file << i << "," 
                 << spatial_coords[i][0] << "," 
                 << spatial_coords[i][1] << "," 
                 << (is_sketched[i] ? 1 : 0) << "\n";
    }

    out_file.close();
    std::cout << "Successfully exported Visium results to " << filename << "\n\n";
}

std::pair<double, double> RunBaselineVisiumPipeline(const std::string& experiment_name,
                       const std::vector<std::vector<int>>& sparse_genes,
                       const std::vector<std::vector<double>>& spatial_coords,
                       int k) {
                     
    std::cout << "--- Running " << experiment_name << " (N=" << spatial_coords.size() << ", k=" << k << ") ---\n";
    
    SpatialSketcher* my_sketcher = new BaselineSketching(sparse_genes, spatial_coords);
    std::vector<std::size_t> final_sketch = my_sketcher->ComputeSketch(k);

    // Metrics
    double rna_hausdorff = ComputeTranscriptomicHausdorff(sparse_genes, final_sketch);
    double spatial_hausdorff = ComputeCoordinateHausdorff(spatial_coords, final_sketch);

    std::cout << "     [Metric] Transcriptomic Hausdorff (95%): " << rna_hausdorff << "\n";
    std::cout << "     [Metric] Coordinate Hausdorff (95%): " << spatial_hausdorff << "\n";
    
    std::string filename = "bin/" + experiment_name + "_results.csv";
    // Call the clean Visium exporter
    ExportVisiumResultsToCSV(filename, spatial_coords, final_sketch);
    
    delete my_sketcher;

    return {rna_hausdorff, spatial_hausdorff};
}

std::pair<double, double> RunQuadtreeVisiumPipeline(const std::string& experiment_name,
                       const std::vector<std::vector<int>>& sparse_genes,
                       const std::vector<std::vector<double>>& spatial_coords,
                       int k, int capacity=-1, int budget=-1) {
                     
    std::cout << "--- Running " << experiment_name << " (N=" << spatial_coords.size() << ", k=" << k << ") ---\n";
    
    SpatialSketcher* my_sketcher = new QuadTreeSketching(sparse_genes, spatial_coords, capacity, budget);
    std::vector<std::size_t> final_sketch = my_sketcher->ComputeSketch(k);

    // EVALUATION BLOCK 
    double rna_hausdorff = ComputeTranscriptomicHausdorff(sparse_genes, final_sketch);
    double spatial_hausdorff = ComputeCoordinateHausdorff(spatial_coords, final_sketch);
    
    std::cout << "     [Metric] Transcriptomic Hausdorff (95%): " << rna_hausdorff << "\n";
    std::cout << "     [Metric] Coordinate Hausdorff (95%): " << spatial_hausdorff << "\n";
    
    std::string filename = "bin/" + experiment_name + "_results.csv";
    // Call the clean Visium exporter
    ExportVisiumResultsToCSV(filename, spatial_coords, final_sketch);
    
    delete my_sketcher;

    return {rna_hausdorff, spatial_hausdorff};
}


int main() {
    
    // // EXPERIMENT 1: 5-Bin Toy Example
    // std::vector<std::vector<int>> toy_genes = {{1, 2}, {1, 2, 3}, {4, 5}, {4}, {1, 5}};
    // std::vector<std::vector<double>> toy_coords = {{0.0, 0.0}, {0.0, 2.0}, {10.0, 10.0}, {10.0, 8.0}, {5.0, 5.0}};
    
    // RunPipeline("toy_dataset", toy_genes, toy_coords, 3);


    // // EXPERIMENT 2: Generalized Grid (400 bins)
    // std::vector<std::vector<int>> grid_genes;
    // std::vector<std::vector<double>> grid_coords;
    
    // std::mt19937 sim_rng(123); 
    // std::uniform_real_distribution<double> prob_dist(0.0, 1.0);
    // std::uniform_int_distribution<int> noise_gene_dist(10, 100); // Background noise genes
    
    // // Create a 20x20 grid of tissue
    // for (int x = 0; x < 20; ++x) {
    //     for (int y = 0; y < 20; ++y) {
    //         grid_coords.push_back({static_cast<double>(x), static_cast<double>(y)});
    //         std::vector<int> current_genes;
            
    //         // Define three hypothetical transcriptomic categories:
    //         // duct core (epithelial), duct boundary (mixed), and outer stroma.
    //         const int dist2 = (x - 10) * (x - 10) + (y - 10) * (y - 10);

    //         if (dist2 < 16) {
    //             // Duct core: mostly epithelial
    //             if (prob_dist(sim_rng) < 0.8) current_genes.push_back(1);
    //             if (prob_dist(sim_rng) < 0.8) current_genes.push_back(2);
    //             if (prob_dist(sim_rng) < 0.8) current_genes.push_back(3);
    //         } else if (dist2 < 36) {
    //             // Boundary ring: mixed expression
    //             if (prob_dist(sim_rng) < 0.7) current_genes.push_back(1);
    //             if (prob_dist(sim_rng) < 0.7) current_genes.push_back(2);
    //             if (prob_dist(sim_rng) < 0.7) current_genes.push_back(4);
    //             if (prob_dist(sim_rng) < 0.7) current_genes.push_back(5);
    //         } else {
    //             // Outer tissue: mostly stromal
    //             if (prob_dist(sim_rng) < 0.8) current_genes.push_back(4);
    //             if (prob_dist(sim_rng) < 0.8) current_genes.push_back(5);
    //             if (prob_dist(sim_rng) < 0.8) current_genes.push_back(6);
    //         }
            
    //         // Add some background noise
    //         current_genes.push_back(noise_gene_dist(sim_rng));
    //         if (prob_dist(sim_rng) < 0.5) {
    //             current_genes.push_back(noise_gene_dist(sim_rng));
    //         }
            
    //         grid_genes.push_back(current_genes);
    //     }
    // }
    
    // RunPipeline("grid_dataset", grid_genes, grid_coords, 100);

    std::cout << "Loading Visium HD Data..." << std::endl;
    
    // Load the preprocessed data
    auto visium_coords = LoadSpatialCoords("data/visium_spatial_coords.csv");
    auto visium_genes = LoadSparseGenes("data/visium_sparse_genes.txt");
    
    if (visium_coords.empty() || visium_genes.empty()) {
        std::cerr << "Failed to load data files. Check paths." << std::endl;
        return 1;
    }

    std::cout << "Loaded " << visium_coords.size() << " spatial bins." << std::endl;

    // Run the specific Visium pipeline!
    // RunQuadtreeVisiumPipeline("visium_quadtree_1000", visium_genes, visium_coords, 1000, 10000, 21);
    // RunBaselineVisiumPipeline("visium_baseline_1000", visium_genes, visium_coords, 1000);

    // New parameter sweep for QuadTree:
    // std::vector<int> capacities = {5000, 10000, 15000, 20000, 30000, 40000, 60000, 80000, 120000, 240000};
    // std::vector<double> multipliers = {1.0, 1.25, 2.0, 3.0, 5.0};
    
    // int k = 1000;
    // double N = static_cast<double>(visium_coords.size()); // ~479863

    // std::cout << "\n=== STARTING DYNAMIC PARAMETER SWEEP ===\n";
    // for (int cap : capacities) {
    //     for (double mult : multipliers) {
            
    //         // Calculate the dynamic budget
    //         double base_budget = (k * static_cast<double>(cap)) / N;
    //         int dynamic_budget = static_cast<int>(std::ceil(base_budget * mult));
    //         if (dynamic_budget < 1) dynamic_budget = 1;

    //         std::string exp_name = "sweep_cap" + std::to_string(cap) + "_mult" + std::to_string(mult);
            
    //         RunQuadtreeVisiumPipeline(exp_name, visium_genes, visium_coords, k, cap, dynamic_budget);
    //     }
    // }

    // Quadtree experiments
    // RunQuadtreeVisiumPipeline("quadtree_500000_1005", visium_genes, visium_coords, 1000, 500000, 1005);
    // RunQuadtreeVisiumPipeline("quadtree_35000_70", visium_genes, visium_coords, 1000, 30000, 70);
    // RunQuadtreeVisiumPipeline("quadtree_30000_60", visium_genes, visium_coords, 1000, 25000, 40);
    // RunQuadtreeVisiumPipeline("quadtree_25000_50", visium_genes, visium_coords, 1000, 15000, 15);
    // RunQuadtreeVisiumPipeline("quadtree_20000_40", visium_genes, visium_coords, 1000, 10000, 10);

    // =================== TODO RYAN ==================
    // Do any paramterization on the voronni graph you feel necessary
    // TODO: Can you run the following parameter sweep to compare baseline, quadtree, and voronoi across $k$ values
    // Can you then take these results and make a plot of coordinate hausdorff and transcriptomic hausdorff across $k$?
    // In the benchmarking paper, they have a 2 table figure where the top is coordinate hausdorff and the bottom is transcriptomic hausdorff 
    // Each figure should have the sampling fraction (k/N) on the x axis and the hausdorff distance on the y axis
    // Then we have 3 lines plots (one for each method) in each figure showing how the methods hausdorff changes across $k$ 
    

    std::vector<int> k_values = {5000, 10000, 25000, 50000, 75000, 100000};
    int optimal_capacity = 10000;
    double N = static_cast<double>(visium_coords.size());

    std::cout << "\n=== STARTING K-SCALING BENCHMARK ===\n";

    // CSVs for visualization
    std::ofstream scaling_file("scaling_metrics.csv");
    scaling_file << "Method,k,Coordinate_Hausdorff,Transcriptomic_Hausdorff\n";

    for (int k : k_values) {
        std::cout << "\nEvaluating k = " << k << "...\n";
        auto baseline_metrics = RunBaselineVisiumPipeline("baseline_k" + std::to_string(k), visium_genes, visium_coords, k);
        scaling_file << "Baseline," << k << "," << baseline_metrics.second << "," << baseline_metrics.first << "\n";

        // Calculate the dynamic budget based on the current k
        double base_budget = (k * static_cast<double>(optimal_capacity)) / N;
        int dynamic_budget = static_cast<int>(std::ceil(base_budget * 1.0)); // 1.0x Multiplier
        if (dynamic_budget < 1) dynamic_budget = 1;

        auto quadtree_metrics = RunQuadtreeVisiumPipeline("quadtree_k" + std::to_string(k), visium_genes, visium_coords, k, optimal_capacity, dynamic_budget);
        scaling_file << "QuadTree," << k << "," << quadtree_metrics.second << "," << quadtree_metrics.first << "\n";

        // Should define a RunVoronoiVisiumPipeline function that matches the signature of RunQuadtreeVisiumPipeline
        auto voronoi_metrics = RunVoronoiVisiumPipeline("voronoi_k" + std::to_string(k), visium_genes, visium_coords, k);
        scaling_file << "Voronoi," << k << "," << voronoi_metrics.second << "," << voronoi_metrics.first << "\n";
        
    }
    scaling_file.close();
    std::cout << "\nScaling benchmark complete! Data saved to scaling_metrics.csv\n";

    return 0;
}