# Epigenome analysis pipeline Part 1


Define regulatory elements (2kb)
- Roadmap DHS Promoter (45963 regions)
- Roadmap DHS Enhancer (619368 regions)
- ATAC peak Bone marrow (311510 regions) : Corces et al.
- ATAT peak Bone marrow and Thymus combined (423488 regions) : unpublished data

Exclusion of the chromosomes
- For the clustering of multiple individuals : exclude chrX, Y, M
- For the clustering of cells within one individual : exclude chrM

Analysis pipeline (workflow management using Snakemake)
- Input: bamfile, remove duplicate and supplementary reads
- Generate count matrix (Deeptool)
- Inference of motif activity (Rscript)
- Supervised cell-type prediction (Rscript)
- Unsupervised clustering (Rscript)
- Trajectory analysis (Rscript)
- Output: cell-type label, clustering plot, motif activity, trajectory plot

Example datasets
- Cell line mixture
- Primary cell data
