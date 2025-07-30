### Analysis of Control and Treated Breast Cancer: Utilization of Differential Gene Expression (DGE)

This project presents a reproducible RNA-Seq pipeline for breast cancer differential gene expression (DGE) analysis, designed with HPC compatibility. It explores transcriptomic data from breast cancer tissue samples to identify genes associated with disease progression, molecular subtypes, and pathway enrichment. Leveraging NextFlow for reproducibility and DESeq2 for statistical rigor, the pipeline integrates preprocessing, normalization, differential testing, and visualization.

##  Objectives
- Detect significantly differentially expressed genes (DEGs) between tumor subtypes and normal tissue.
- Visualize expression profiles to uncover biological signals.
- Perform enrichment analysis to infer functional pathways and disease mechanisms.
- Demonstrate reproducible, scalable bioinformatics using NextFlow

## Files in this Repository

| File | Description |
| `Report` | Final report including the project summary, methods, results, discussion, and code |
| `PCA_plot` | PCA Plot
| `ma_plot_preshrinkage.png` | MA Plot (Pre-Shrinkage)
| `ma_plot_postshrinkage.png` | MA Plot (Post-Shrinkage)
| `dispersion_plot.png` | Dispersion-by-mean Plot
| `pvalue_histogram.png` | Raw P-value Histogram
| `volcano_plot.png` | Volcano Scatterplot
| `sample_mapping_metrics` | Sample Mapping Metrics
| `dge_top10` | Top 10 Significant DGEs
| `dge_summary_counts` | DGE Summary (Upregulated vs Downregulated Counts)
| `deseq2_results_full` | Full DESeq2 Output
| `workflow_execution_screenshots` | Screenshot of Nextflow Execution

## Key Methods

# Sample Preparation & Data Acquisition
-Samples: RNA-seq data from control and NRDE2 RNAi-treated breast cancer cell lines.
-Sequencing Type: Single-end reads from Illumina NextSeq platform.
-Metadata: Categorized into control and treated groups (3 replicates each).
-Reference Files: Human reference transcriptome (FASTA) and GTF annotation from ENSEMBL.

# Data Processing Pipeline
-Workflow Tool: nf-core/rnaseq pipeline (v3.14.0) executed via Nextflow.

-HPC Execution: Submitted as SLURM job with custom config (rna.json) and sample sheet.

-Read Trimming & QC:
---TrimGalore! and FastQC for adapter and quality trimming.
---Aggregated using MultiQC for overall quality assessment.

-Transcript Quantification:
---Performed using Salmon, which provides quantification via quant.sf files.

# Gene-Level Analysis
-Transcript-to-Gene Mapping:
---Used tximport to summarize transcript-level counts into gene-level estimates.

-Preprocessing in R:
---Filtering out low-expression genes (<10 total reads).
---Normalization via DESeq2 using the median-of-ratios method.
---Variance-stabilizing transformation with rlog for PCA.

# Differential Expression Testing
-Statistical Model:
---Negative binomial model via DESeq2.
---Log2 Fold Change (LFC) estimation enhanced with apeglm shrinkage.
---Adjusted p-values computed using Benjamini-Hochberg correction (FDR < 0.05).

# Visualization & Interpretation
-Plots Generated:
---PCA plot (using plotPCA) to visualize group separation.
---MA plots (pre- and post-shrinkage) to assess expression change distribution.
---Volcano plot to highlight significant DEGs.
---Dispersion-by-mean plot and raw p-value histogram for QC and statistical diagnostics.

# Key Output Tables:
-Mapping metrics per sample.

-Top significant DEGs.

-Full DESeq2 results table.

## Results

| Gene Symbol | log2 Fold Change | Adjusted p-value | Regulation   |
|-------------|------------------|------------------|--------------|
| DDX58       | +2.67            | 3.5e-06          | Upregulated  |
| TP53INP1    | +2.03            | 1.2e-04          | Upregulated  |
| ZNF703      | –2.94            | 4.9e-05          | Downregulated|
| BCL11A      | –2.32            | 6.7e-04          | Downregulated|
| CCL5        | +1.87            | 2.9e-03          | Upregulated  |






