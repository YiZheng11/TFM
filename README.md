# Mitochondrial Sex Dimorphism (TFM Code)
This repository contains the shell scripts, R scripts, and Jupyter notebooks used in my Master’s thesis to analyse sexual dimorphism in mitochondrial gene expression in mouse liver and human liver/heart. The code implements a reproducible pipeline from read preprocessing and counting (mouse RNA‑seq) to differential expression and functional analyses on MitoCarta‑annotated genes in mouse and GTEx datasets.

Raw sequencing data (mouse) and GTEx files are not included for privacy and storage reasons. The repository focuses on analysis scripts and how to regenerate the main summary tables and figures from the corresponding input files.

## Repository structure
```
.
├── 00_read_distribution.R          # QC summaries and read distribution plots for count data
├── 01_deseq.R                      # Core DESeq2 differential expression for mouse liver
├── 01A_tsv2excel.ipynb             # Export of DE/summary tables to Excel
├── 01B_GO_classification.R         # GO classification/enrichment of significant genes
├── 02_log2FC.R                     # Log2FC extraction and visualisation for selected gene sets
├── 03_GSEA.R                       # GSEA on ranked gene lists (mouse liver)
├── 04A_gtex_mitocarta.ipynb        # GTEx: import, filtering and MitoCarta subsetting
├── 04B_gtex_distribution.R         # GTEx: sample distribution and exploratory plots
├── 04C_deseq_gtex.R                # DESeq2 analyses for human liver and heart
├── 05A_comparative_analysis_genes.R      # Cross‑species comparison at gene level
├── 05B_comparative_analysis_function.R  # Cross‑species comparison at pathway level
├── Mm_gtf_filt.sh                  # Mouse GTF preprocessing restricted to MitoCarta genes
├── htseq_gtf_filt.sh               # HTSeq counting script using filtered GTF
└── trim.sh                         # Read trimming and basic QC for mouse FASTQ files
```
Directories for data, results and figures should be created locally; adjust paths inside the scripts to match your environment.

## Main analyses

**Mouse liver (SKA111/SKA113)**
- Preprocessing on the CBMSO cluster: `trim.sh, Mm_gtf_filt.sh, htseq_gtf_filt.sh`
- QC and differential expression: `00_read_distribution.R, 01_deseq.R`
- Post‑processing and functional analyses:
  - `01A_tsv2excel.ipynb` - export DESeq2 tables to Excel.
  - `01B_GO_classification.R` - GO enrichment plots.
  - `02_log2FC.R` - log2FC extraction and visualisation for selected gene sets.
  - `03_GSEA.R` - GSEA on ranked gene lists.

**Human liver and heart (GTEx)**
- `04A_gtex_mitocarta.ipynb` – import GTEx expression/metadata and subset to MitoCarta genes.
- `04B_gtex_distribution.R` – age/sex distributions and exploratory plots.
- `04C_deseq_gtex.R` – sex- and age-stratified DESeq2 analyses for liver and heart (atrial appendage and left ventricle).
- `03_GSEA.R` - GSEA on ranked gene lists.

**Comparative mouse–human liver analyses**
- `05A_comparative_analysis_genes.R` – orthologous mitochondrial DE genes and overlaps between mouse and human.
- `05B_comparative_analysis_function.R` – shared and species-specific mitochondrial/metabolic pathways at the functional level.

## Usage
1. Mouse RNA‑seq preprocessing

    Run `trim.sh, Mm_gtf_filt.sh and htseq_gtf_filt.sh` on your cluster, adjusting paths to FASTQ, BAM and GTF files. These scripts generate the gene‑level count tables used as input for DESeq2.

2. Mouse liver DE and functional analysis

- `00_read_distribution.R` to summarise counts and produce QC plots.
- `01_deseq.R` to run DESeq2 and obtain normalised counts and DEGs.
- `01A_tsv2excel.ipynb, 01B_GO_classification.R, 02_log2FC.R and 03_GSEA.R` to generate the final tables and figures used in the thesis.

3. GTEx liver and heart analyses

- `04A_gtex_mitocarta.ipynb` to load GTEx, filter mitochondrial genes and prepare count matrices.
- `04B_gtex_distribution.R` for basic exploratory analyses by sex and age.
- `04C_deseq_gtex.R` for sex/age differential expression in human liver and heart.

4. Cross‑species comparison

- `05A_comparative_analysis_genes.R` and `05B_comparative_analysis_function.R` to reproduce the gene‑ and pathway‑level mouse–human comparisons presented in the TFM.

 ## Reproducibility
Full reproduction of the analyses requires access to:
- Mouse FASTQ/BAM files and the corresponding sample metadata.
- GTEx expression data and metadata for liver, heart atrial appendage and left ventricle.
- Human and mouse MitoCarta 3.0 gene lists.

Given equivalent inputs and correct paths, running the scripts in the order described above will regenerate the main count tables, DESeq2 results, GSEA outputs and most figures reported in the thesis.
