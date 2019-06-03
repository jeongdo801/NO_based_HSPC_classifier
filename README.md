# Epigenome analysis pipeline Part 1


Define regulatory elements (2kb)
- Roadmap DHS Enhancer (619368 regions)
- ATAC peak Bone marrow (311510 regions) : Corces et al.
- ATAT peak Bone marrow and Thymus combined (423488 regions) : unpublished data


Exclusion of the chromosomes
- For the analysis of multiple individuals : exclude chrX, Y, M (default)
- For the analysis of cells within one individual : exclude chrY, M


Analysis pipeline (workflow management using Snakemake)
- Input: bamfile 
- Preprocessing to remove duplicate and supplementary reads (samtools, biobambam)
- Generate count matrix (Deeptool)
- Inference of motif activity (Rscript)
- Supervised cell-type prediction (Rscript)
- Output: motif activity per single-cell, predicted cell-type label


Example datasets
- Cell line data (LCL, K562, RPE1)
- Primary cell data


Further ideas
- Unsupervised clustering (Rscript)
- Trajectory analysis (Rscript)
