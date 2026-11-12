# Software

Versions are those reported in the Methods of the manuscript (Sections 4.2 and 5.2). "n/r" = version not reported.

## Phylogenomics and genomic similarity (Section 4.2)

| Step | Software | Version | Folder |
| --- | --- | --- | --- |
| Genome download | NCBI Datasets (`datasets`) | n/r | `01_download_genomes` |
| Completeness | BUSCO (`aspergillus_odb12`, Miniprot) | 5.8.2 | `02_busco_genomes` |
| Orthogroups | OrthoFinder | 2.5.5 | `03_phylogenomics` |
| Similarity search | DIAMOND | 2.2.2 | `03_phylogenomics` |
| Alignment | MAFFT | 7.526 | `03_phylogenomics` |
| ML tree / model selection | IQ-TREE / ModelFinder | 3.1.2 | `03_phylogenomics` |
| Species delimitation (set-up) | BEAST2 + STACEY add-on, jModelTest | n/r | `03_phylogenomics` |
| Tree display | FigTree | 1.4.4 | - |
| Alignment-free ANI | Mash (k = 21, sketch = 10,000) | 2.3 | `04_genomic_similarity` |
| Alignment-based ANI | pyANI, method ANIm (MUMmer 3.23) | 0.2.13.1 | `04_genomic_similarity` |
| Heatmaps | R, pheatmap, ggplot2 | n/r | - |
| Assembly purity | seqkit, BLAST+ (nt) | n/r | `05_assembly_purity` |

## Transcriptome (Section 5.2)

| Step | Software | Version | Folder |
| --- | --- | --- | --- |
| Read QC | FastQC | 0.12.1 | - |
| Adapter/quality trimming, k-mer filtering | BBTools (BBDuk) | 38.96 | `01_read_preprocessing` |
| rRNA removal | SortMeRNA (database smr v4.3) | 8.0.1 | `01_read_preprocessing` |
| Assembly | Trinity | 2.15.2 | `02_assembly_and_orfs` |
| Redundancy / overlap | CD-HIT (`cd-hit-est`, `cd-hit-est-2d`) | n/r | `02_assembly_and_orfs` |
| Completeness | BUSCO (`viridiplantae_odb12`, `liliopsida_odb12`) | 6.1.0 | `02_assembly_and_orfs` |
| ORF prediction | TransDecoder | 7.0.0 | `02_assembly_and_orfs` |
| Coding potential | CPC2 | n/r | `02_assembly_and_orfs` |
| Quantification | Salmon | 1.8.0 | `03_quantification` |
| Import of counts | tximport (R) | n/r | `04_differential_expression` |
| Differential expression | DESeq2 (R), `filterByExpr` from edgeR, `lfcShrink` with apeglm | DESeq2 1.44.0; edgeR 3.19 | `04_differential_expression` |
| Cross-check | limma-voom | n/r | `04_differential_expression` |
| Exploratory plots | pheatmap, viridis, ggplot2 | n/r | `04_differential_expression` |
| Domain annotation / GO | InterProScan | 5.76 | `05_annotation_and_enrichment` |
| Enrichment | clusterProfiler, GO.db (R) | n/r | `05_annotation_and_enrichment` |
| Variant calling | samtools, bcftools, FreeBayes | n/r | `06_variant_calling` |
| Validation | BLAST+, seqtk, bowtie2, samtools, IGV | n/r | `07_validation_vs_asparagus` |
| PPI network | STRING, NetworkX, Cytoscape | 12.0; 3.7; 3.10.5 | - |
