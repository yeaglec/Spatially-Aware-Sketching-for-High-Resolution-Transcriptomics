#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <tuple>

#include "model.cpp"

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

void ExportResultsToCSV(const std::string& filename, 
                        const std::vector<std::vector<int>>& sparse_genes,
                        const std::vector<std::vector<double>>& spatial_coords,
                        const std::vector<std::size_t>& sketch_indices) {
    
    std::ofstream out_file(filename);
    if (!out_file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << " for writing.\n";
        return;
    }

    // Create a boolean mask to easily tag which points were sketched
    std::vector<bool> is_sketched(spatial_coords.size(), false);
    for (std::size_t idx : sketch_indices) {
        is_sketched[idx] = true;
    }

    // Write header
    out_file << "Original_Index,X_Coordinate,Y_Coordinate,Epithelial_Markers,Stromal_Markers,Dominant_Program,Is_Sketched\n";

    // Write all data
    for (std::size_t i = 0; i < spatial_coords.size(); ++i) {
        const auto [epi_count, stroma_count, dominant_program] = ExportToyMarkers(sparse_genes[i]);
        out_file << i << "," 
                 << spatial_coords[i][0] << "," 
                 << spatial_coords[i][1] << "," 
                 << epi_count << ","
                 << stroma_count << ","
                 << dominant_program << ","
                 << (is_sketched[i] ? 1 : 0) << "\n";
    }

    out_file.close();
    std::cout << "Successfully exported results to " << filename << "\n\n";
}

// PIPELINE RUNNER
void RunPipeline(const std::string& experiment_name,
                 const std::vector<std::vector<int>>& sparse_genes,
                 const std::vector<std::vector<double>>& spatial_coords,
                 int k) {
                     
    std::cout << "--- Running " << experiment_name << " (N=" << spatial_coords.size() << ", k=" << k << ") ---\n";
    
    SpatialSketcher* my_sketcher = new BaselineSketching(sparse_genes, spatial_coords);
    std::vector<std::size_t> final_sketch = my_sketcher->ComputeSketch(k);
    
    std::string filename = experiment_name + "_results.csv";
    ExportResultsToCSV(filename, sparse_genes, spatial_coords, final_sketch);
    
    delete my_sketcher;
}


int main() {
    
    // EXPERIMENT 1: 5-Bin Toy Example
    std::vector<std::vector<int>> toy_genes = {{1, 2}, {1, 2, 3}, {4, 5}, {4}, {1, 5}};
    std::vector<std::vector<double>> toy_coords = {{0.0, 0.0}, {0.0, 2.0}, {10.0, 10.0}, {10.0, 8.0}, {5.0, 5.0}};
    
    RunPipeline("toy_dataset", toy_genes, toy_coords, 3);


    // EXPERIMENT 2: Generalized Grid (400 bins)
    std::vector<std::vector<int>> grid_genes;
    std::vector<std::vector<double>> grid_coords;
    
    std::mt19937 sim_rng(123); 
    std::uniform_real_distribution<double> prob_dist(0.0, 1.0);
    std::uniform_int_distribution<int> noise_gene_dist(10, 100); // Background noise genes
    
    // Create a 20x20 grid of tissue
    for (int x = 0; x < 20; ++x) {
        for (int y = 0; y < 20; ++y) {
            grid_coords.push_back({static_cast<double>(x), static_cast<double>(y)});
            std::vector<int> current_genes;
            
            // Define three hypothetical transcriptomic categories:
            // duct core (epithelial), duct boundary (mixed), and outer stroma.
            const int dist2 = (x - 10) * (x - 10) + (y - 10) * (y - 10);

            if (dist2 < 16) {
                // Duct core: mostly epithelial
                if (prob_dist(sim_rng) < 0.8) current_genes.push_back(1);
                if (prob_dist(sim_rng) < 0.8) current_genes.push_back(2);
                if (prob_dist(sim_rng) < 0.8) current_genes.push_back(3);
            } else if (dist2 < 36) {
                // Boundary ring: mixed expression
                if (prob_dist(sim_rng) < 0.7) current_genes.push_back(1);
                if (prob_dist(sim_rng) < 0.7) current_genes.push_back(2);
                if (prob_dist(sim_rng) < 0.7) current_genes.push_back(4);
                if (prob_dist(sim_rng) < 0.7) current_genes.push_back(5);
            } else {
                // Outer tissue: mostly stromal
                if (prob_dist(sim_rng) < 0.8) current_genes.push_back(4);
                if (prob_dist(sim_rng) < 0.8) current_genes.push_back(5);
                if (prob_dist(sim_rng) < 0.8) current_genes.push_back(6);
            }
            
            // Add 1 or 2 random background noise genes to EVERY bin to ensure uniqueness
            current_genes.push_back(noise_gene_dist(sim_rng));
            if (prob_dist(sim_rng) < 0.5) {
                current_genes.push_back(noise_gene_dist(sim_rng));
            }
            
            grid_genes.push_back(current_genes);
        }
    }
    
    
    // We have 400 bins. Let's see what happens if we only sketch 15 of them.
    RunPipeline("grid_dataset", grid_genes, grid_coords, 100);

    return 0;
}