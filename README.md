# Seurat-scRNA-seq-Full-Analysis-PBMC-3k-
Automated End to End Single Cell RNA Seq Analysis Pipeline Quality Control Unsupervised Clustering and Interactive 3D Visualization PBMC 3k
Automated Single-Cell RNA-Seq Analysis Pipeline (PBMC 3k)
🧬 Overview

This repository contains an end-to-end R pipeline for the analysis of the Peripheral Blood Mononuclear Cell (PBMC) 3k dataset from 10x Genomics.

The script automates the entire process, from raw data acquisition and Quality Control (QC) through advanced dimensionality reduction, unsupervised clustering, marker gene detection, and final cell-type annotation. A key feature is the generation of interactive 2D and 3D visualization files for dynamic data exploration.
Key Features

    Automated Setup: Checks and installs necessary R packages (Seurat, plotly).

    Data Source: Automatically downloads the raw data matrix from 10x Genomics.

    Processing: Includes best-practice steps: QC, Normalization, Scaling, PCA, UMAP.

    Interactive Results: Generates zoomable, hoverable HTML plots, including a 3D UMAP cube.

    Annotation: Automatically renames clusters (0, 1, 2...) to their definitive biological names (T-Cells, B-Cells, Monocytes, etc.).

    Output Management: Saves all results (plots, data, final R object) directly to the execution directory.

🚀 Getting Started
Prerequisites

    R and RStudio: Ensure you have R (version 4.0 or higher) and RStudio installed.

    Internet Access: Required to download the raw data and install packages.

    Pandoc: Required by the htmlwidgets package to save the interactive HTML files. (If using RStudio, this is usually included.)

Installation and Execution

The script is designed to be run from the RStudio Console or as a single job.

    Download the Script: Clone this repository and save the R script (Automated End to End Single Cell RNA Seq Analysis Pipeline...R) to your local machine.

    Run the Script: Open the script in RStudio and press Ctrl+A (select all) followed by Ctrl+Enter (run selected lines).

The script will automatically handle package installation, data download, and set the working directory to the folder where the script is run.
📂 Results and Outputs

The script saves the following files to your execution directory:

Filename,Type,Description
09_Final_Annotated_UMAP.png,PNG,"Final, publication-ready figure with cell types clearly labeled (e.g., Naive CD4 T)."
06_Interactive_3D_UMAP.html,HTML,Interactive 3D plot of the clusters. Open in any web browser to rotate and explore the data.
07_Cluster_Markers.csv,CSV,Table listing the top differentially expressed genes (marker genes) for every cluster.
04_Interactive_UMAP.html,HTML,Standard 2D interactive UMAP plot (hover to view cell barcodes/clusters).
pbmc3k_annotated_final.rds,RDS,"The final saved R object (Seurat object) containing all the processed data, annotations, and UMAP coordinates."
🔬 Scientific Interpretation

The analysis successfully identified 9 distinct cell populations from the PBMC 3k dataset, consistent with published single-cell literature.

The final annotations applied in Step 8 of the script are:
Cluster ID,Annotated Cell Type,Defining Marker Genes
"0, 1",Naive/Memory CD4+ T,"IL7R, CCR7, S100A4"
2,CD14+ Monocytes,"CD14, LYZ"
3,B Cells,"MS4A1, CD79A"
4,CD8+ T Cells,CD8A
6,NK Cells,"GNLY, NKG7"
"7, 8",DC / Platelets,"FCER1A, PPBP"
Data Exploration

To further explore the results, load the
library(Seurat)
pbmc <- readRDS("pbmc3k_annotated_final.rds")

# Generate a final DotPlot summarizing key markers
features_to_plot <- c("CD14", "LYZ", "MS4A1", "CD3E", "CD8A", "GNLY", "PPBP")
DotPlot(pbmc, features = features_to_plot) + RotatedAxis()🛠️ Dependencies

This script requires the following R packages, which it attempts to install automatically:

    Seurat (Core scRNA-seq analysis)

    dplyr (Data manipulation)

    ggplot2 (Static plotting)

    patchwork (Plot arrangement)

    plotly (Interactive visualization)

    htmlwidgets (Saving interactive output)
