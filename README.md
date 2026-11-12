# Resources used for the article *Aspergillus welwitschiae* and a New Cryptic Species

Scripts and supplementary data accompanying the manuscript **"Untangling the Nigri Section: A Phylogenomic Reassessment of
*Aspergillus welwitschiae* and a New Cryptic Species"** and the companion transcriptomic analysis of *Agave sisalana* (sisal) infected by
*A. welwitschiae* (red stem rot), prepared for the MSc qualification at the Universidade Federal de Minas Gerais (UFMG).

## Repository layout

| Folder | Content | Manuscript section |
| --- | --- | --- |
| `01_aspergillus_genomics/` | Genome download, BUSCO quality control, OrthoFinder + IQ-TREE phylogenomics, Mash and pyANI (ANIm) similarity, assembly purity check, marker search | 4.2 |
| `02_agave_transcriptome/` | Read pre-processing, de novo assembly, quantification, differential expression, exploratory analyses, GO enrichment, variant calling, validation | 5.2 |
| `supplementary_tables/` | Large supplementary tables that are too big for the manuscript (Tables S3, S4, ...) | Supplementary Material |
| `legacy_preliminary_analyses/` | Earlier approaches that are not part of the reported results | - |
| `docs/` | `SOFTWARE.md` (tools and versions) and `PATHS.md` (directory placeholders) | - |

Each analysis step has its own numbered folder; scripts are meant to be read and run in numerical order.

## Directory placeholders

Paths that depend on the local environment were replaced by placeholders in square brackets, for example `[diretório_transcriptoma]`.
Replace each placeholder with the corresponding path on your machine before running a script. The full legend is in
[`docs/PATHS.md`](docs/PATHS.md); to list every placeholder used in a script:

```bash
grep -o "\[diretório_[a-z_0-9]*\]" -r . | sort | uniq -c
```

## Software

Tools and versions are listed in [`docs/SOFTWARE.md`](docs/SOFTWARE.md).

## Supplementary tables stored here

See [`supplementary_tables/README.md`](supplementary_tables/README.md).
