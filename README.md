# Spatially-Aware-Sketching-for-High-Resolution-Transcriptomics
Standard transcriptomic sketching algorithms (MinHash) are effective at downsampling in high-dimensional gene-expression space, but are unable to handle high-resolution spatial transcriptomics (VisiumHD). When applied to these datasets, standard sketching algos induce geometric distortion in tissue structures.

## Dataset

We use the 10x Genomics Visium HD 3' Gene Expression Library, Human Pancreatic Cancer (Fresh Frozen).
The raw data is available at [https://www.10xgenomics.com/datasets/visium-hd-three-prime-human-pancreatic-cancer-fresh-frozen]

Download the **Binned outputs (all bin levels)** tar ball, uncompress it, then move the files into the **data** directory. Optional: You can download the **Spatial** directory as well for visualization purposes.

Assuming data directory is set up as intended, run **pancreatic_data_processing.ipynb** for data pre-processing. This script creates **visium_sparse_genes.txt** and **visium_spatial_coords.csv** in the data directory. These are the gene expression matrix (sparse) and spatial coordinates matrix used in the sketch algorithms. 

For eval, run the pipeline with the function calls ComputeTranscriptomicHausdorff(sparse_genes, final_sketch) and ComputeCoordinateHausdorff(spatial_coords, final_sketch) in main.cpp right before you export the results to CSV.


