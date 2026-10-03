# U1SiteScore

Code and processed inputs for analyzing the genomic distribution of U1 snRNP
binding-site scores (U1SS) and their association with human gene expression.
The scripts in this repository accompany the associated journal publication.

## Repository contents

| File | Description |
| --- | --- |
| `map_5ss_to_genome.py` | Maps 9-mer 5′ splice-site scores across a reference genome and writes strand-specific bedGraph tracks. |
| `Find_top-n_nearest_U1sites.R` | Identifies the nearest strand-matched U1SS intervals for each protein-coding gene transcription start site. |
| `Generalized-Additive-Model_for-U1SS.Rmd` | Fits the negative-binomial generalized additive model relating U1SS features to mean expression across ten human tissues. |
| `U1site-distribution-whole-genome.Rmd` | Calculates U1SS density in promoter, exon, intron, and distal intergenic regions. |
| `nearest_U1sitelist.list` | Processed R data file containing the ten nearest strand-matched U1SS intervals for each gene. |
| `u1_summary_down500_sense_maxentscore_higherthan0.txt` | Per-gene summary of downstream sense-oriented U1 sites within 500 bp with a MaxEnt score greater than zero. |

The `.list` file is an R serialized-data file created with `save()`; it is not a
plain-text list. It contains an object named `nearest_U1sitelist`.

## Software requirements

The R analyses were validated with R 4.3 or later. The genome-mapping workflow
requires Python 3 and Biopython.

Install the required R packages:

```r
install.packages(c("dplyr", "ggplot2", "knitr", "mgcv", "rmarkdown"))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c(
  "GenomeInfoDb",
  "GenomicRanges",
  "IRanges",
  "rtracklayer",
  "S4Vectors"
))
```

Install the Python dependency:

```bash
python -m pip install biopython
```

Rendering the R Markdown documents also requires Pandoc. It is included with
RStudio and can also be installed separately.

## Input data

The GAM uses the two processed U1SS files included in this repository and the
following additional file:

| Input | Required columns or contents |
| --- | --- |
| `rna_tissue_consensus.tsv` | Human Protein Atlas tissue-expression table containing `Gene`, `Gene name`, `Tissue`, and `nTPM`. |

Place `rna_tissue_consensus.tsv` in the repository root. Alternatively, set
`U1SS_DATA_DIR` to the directory containing all three GAM input files.

Recreating all processed U1SS inputs from their upstream sources additionally
requires:

| Input | Used by |
| --- | --- |
| Reference-genome FASTA | `map_5ss_to_genome.py` |
| Two-column table of 9-mer sequences and scores | `map_5ss_to_genome.py` |
| Plus- and minus-strand U1SS bigWig tracks | Both R-based genomic analyses |
| Protein-coding gene TSS table with `chr`, `TSS`, `strand`, and `gene_name` | `Find_top-n_nearest_U1sites.R` |
| GENCODE exon/transcript table described below | `U1site-distribution-whole-genome.Rmd` |

The GENCODE table must contain `Gene.ID.version`, `gene_name`, `chr`,
`exon_start`, `exon_end`, `Exon_rank`, `transcript_start`, `transcript_end`, and
`strand`.

Large reference files and genome-browser tracks are not included in this
repository. The R scripts accept alternate locations through environment
variables and report missing inputs before beginning an analysis.

## Analysis workflow

### 1. Map 9-mer scores across the genome

This optional upstream step creates plus- and minus-strand bedGraph files:

```bash
python map_5ss_to_genome.py \
  --fasta /path/to/reference.fa \
  --scores /path/to/u1ss_9mer_scores.tsv \
  --output-prefix results/map_5ss_hg38 \
  --workers 8 \
  --window-size 100 \
  --step-size 10
```

The outputs are `results/map_5ss_hg38_plus.bg` and
`results/map_5ss_hg38_minus.bg`. Conversion to indexed bigWig tracks is an
external preprocessing step and is not performed by this script.

### 2. Find the nearest U1SS intervals

Place the strand-specific bigWig files and gene-coordinate table in one input
directory, then run:

```bash
U1SS_DATA_DIR=/path/to/input \
U1SS_RESULTS_DIR=results \
U1SS_TOP_N=10 \
Rscript Find_top-n_nearest_U1sites.R
```

The default output is `results/nearest_U1sitelist.list`. Individual paths can
instead be set with `U1SS_PLUS_BIGWIG`, `U1SS_MINUS_BIGWIG`,
`U1SS_GENE_TABLE`, and `U1SS_NEAREST_OUTPUT`.

### 3. Fit the generalized additive model

With `nearest_U1sitelist.list`,
`u1_summary_down500_sense_maxentscore_higherthan0.txt`, and
`rna_tissue_consensus.tsv` in the repository root, render:

```bash
Rscript -e 'rmarkdown::render("Generalized-Additive-Model_for-U1SS.Rmd")'
```

To keep the inputs elsewhere:

```bash
U1SS_DATA_DIR=/path/to/input \
U1SS_RESULTS_DIR=results \
Rscript -e 'rmarkdown::render("Generalized-Additive-Model_for-U1SS.Rmd")'
```

The model uses the mean expression across cerebral cortex, heart muscle,
liver, kidney, lung, spleen, skin, adipose tissue, ovary, and small intestine.
It retains genes at or below the 95th expression percentile, samples 10,000
genes for training with random seed 321, and assigns all genes to ten groups by
predicted expression.

Outputs written to `results/`:

- `gam_gene_predictions.tsv`
- `gam_model.rds`
- `gam_model_summary.txt`
- `gam_observed_expression_by_prediction_decile.pdf`

### 4. Calculate genome-wide regional density

Place the GENCODE table and strand-specific bigWig files in one directory and
run:

```bash
U1SS_DATA_DIR=/path/to/input \
U1SS_RESULTS_DIR=results \
Rscript -e 'rmarkdown::render("U1site-distribution-whole-genome.Rmd")'
```


This analysis uses chromosomes 1–22, X, and Y; transcripts with at least five
annotated exons; five non-overlapping 1-kb promoter bins; exon and intron ranks
1–4; later exons and introns; and strand-matched distal intergenic sequence.


## Citation

If you use this repository, please cite [Kim et al., 2026]. DOI to be added.
