# CN Thermoregulatory Network

R scripts used for the analysis of phosphoproteome, TurboID, RNA-seq and Ribo-seq data presented in the Yadav and Heitman manuscript investigating the thermoregulatory network of *Cryptococcus neoformans*.

## Overview

This repository contains R scripts used for analysis and visualization of phosphoproteomics, TurboID proximity-labeling, RNA sequencing, and Ribo-seq datasets generated as part of this study. The scripts include data processing, statistical analysis, differential analysis, and visualization of the datasets. 

---

## Repository Contents

### Phosphoproteomics

The following scripts were used for analysis of phosphoproteomic datasets to identify shared substrates between calcineurin and Yak1, as well as to analyze the impact on phosphorylation due to various truncation alleles of catalytic subunit of calcineurin, Cna1. The scripts were used for differential phosphoproteomic analysis, generating PCA plots, heatmaps, and UpSet plots. The input files are available in the repository "Data" folder.

| Script                                                                                 | Description                                                                                                                  | Input file                                                                                         |
| -------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------|
| DEP2-script-CN-Yak1-phospho-control-analysis-PCA-plot-with-labels.R                    | Generates PCA plots for the calcineurin, yak1 mutant phosphoproteomics dataset.                                              | 10514_SupplementalData_H99_081524-normalized-WT-cna1-yak1-renamed-dup-split.csv                    |
| DEP2-script-CN-Yak1-phospho-manual-groups-analysis.R                                   | Performs analysis of the phosphoproteomics dataset using manually defined experimental groups.                               | 10514_SupplementalData_H99_081524-normalized-WT-cna1-yak1-renamed-dup-split.csv                    |
| DEP2-script-CN-Yak1-phospho-manual-groups-analysis-heatmap-TurboID-protein-specific.R  | Generates protein-specific phosphoproteomic heatmaps combined with the TurboID analysis.                                     | 10514_SupplementalData_H99_081524-normalized-WT-cna1-yak1-renamed-dup-split.csv                    |
| DEP2-script-CN-truncations_phospho-control-analysis-UpSet-plot.R                       | Generates UpSet plots to visualize overlap among phosphoproteomic datasets from truncation experiments at the protein level. | 10514_SupplementalData_H99_081524-normalized-all-dataset-except-Cnb1-renamed-row-dup-col-split.csv |
| DEP2-script-CN-truncations_phospho-control-analysis-downregulated_Suburst_UpSet-plot.R | Generates SunBurst plots for analysis of hyperphosphorylated or hypophosphorylated proteins at the peptide level.            | 10514_SupplementalData_H99_081524-normalized-all-dataset-except-Cnb1-renamed-row-dup-col-split.csv |

### TurboID Proximity Labeling

The following scripts were used for the analysis and visualization of TurboID proximity-labeling experiments. The results generated from these were also combined with the phosphoproteome data in one of the scripts described above. The input files are available in the repository "Data" folder.

| Script                                                            | Description                                                                                                                                       | Input file                                                             |
| ------------------------------------------------------------------| ------------------------------------------------------------------------------------------------------------------------------------------------- | -----------------------------------------------------------------------|
| DEP2-script-TurboID-data-analysis.R                               | Performs analysis of TurboID proximity-labeling data for initial comparisons.                                                                     | 5675_SupplementalData_112921_normalized-data-TurboID-renamed.csv       |
| DEP2-script-TurboID-data-control_analysis-PCA_heatmap.R           | Performs basic analysis of TurboID data using control method, and visualizes the data including PCA and heatmap visualization.                    | 5675_SupplementalData_112921_normalized-data-TurboID-renamed.csv       |
| DEP2-script-TurboID-data-analysis-Venn-diagram-group-analysis.R   | Provides detailed analysis of various groups and differential enrichment comparing groups at different temperatures and FK506 treated conditions. | 5675_SupplementalData_112921_normalized-data-TurboID-renamed.csv       |

### Ribo-seq

The following scripts were used for the analysis of the Ribo-seq and RNA-seq data using deltaTE tool and identify genes with differential transcription and translation in the calcineurin mutant at 37C. The RFP analysis was done using a custom script and does not use deltaTE tool.

| Script                     | Description                                                        | Input files                                               |
| -------------------------- | ------------------------------------------------------------------ |-----------------------------------------------------------|
| Ribo-seq_RFP_analysis.R    | Analyzes ribosome footprinting (RFP) length from the bam files.    | All bam files are in the "Data/Ribo-seq/bam_files" folder |
| deltaTE_script.R           | Calculates and analyzes changes in translational efficiency (ΔTE). | All files are in the "Data/Ribo-seq" folder               |
| deltaTE_script_PCA_plots.R | Generates PCA plots for translational-efficiency datasets.         | All files are in the "Data/Ribo-seq" folder               |

---

## Software and R Packages

The analyses were performed in **R** using packages including:

* [DEP2](https://bioconductor.org/packages/DEP2/)
* `ggplot2`
* `dplyr`
* `tidyr`
* `readr`
* `pheatmap`
* `ComplexHeatmap`
* `ggrepel`
* `UpSetR`
* `deltaTE`
* Additional packages as specified within individual scripts

Because package requirements may differ between analyses, users should inspect the beginning of each script for the complete list of required packages.

---

## Data

The scripts in this repository were developed for the datasets described in the associated manuscript. The input data files for the described scripts are included in this repository. If additional information or datasets are needed, the users are requested to contact Vikas Yadav directly.

---

## Reproducibility

The scripts were used to generate analyses and figures presented in the associated manuscript.

Most scripts contain paths specific to the computational environment in which the analyses were originally performed. Therefore, users will need to modify file paths and, where necessary, input-file names to match their local environment.

For reproducibility:

1. Install the required R packages.
2. Obtain the appropriate input dataset.
3. Update the input and output paths in the relevant R script.
4. Run the scripts corresponding to the desired analysis.

The scripts are provided to document the computational analyses and facilitate reproduction and further exploration of the results.

---

## Associated Manuscript

Thermoregulation network governing virulence of a critical human fungal pathogen
Vikas Yadav, Joseph Heitman
bioRxiv 2025.08.18.670910; doi: https://doi.org/10.1101/2025.08.18.670910

The manuscript and associated datasets should be consulted for detailed descriptions of the experimental design, biological interpretation, and statistical methods.

---

## Contact

**Vikas Yadav**
[Department of Microbiology and Immunology / University of Minnesota]
yadav@umn.edu

For questions regarding the scripts or analyses, please contact the authors.
