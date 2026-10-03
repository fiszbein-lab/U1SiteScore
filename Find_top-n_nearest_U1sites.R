#!/usr/bin/env Rscript

# Identify the nearest strand-matched U1 splice-site score intervals for each
# protein-coding gene transcription start site (TSS).
#
# Input and output locations default to project-relative directories. They can
# be overridden without editing this script by setting these environment
# variables:
#
#   U1SS_DATA_DIR          Directory containing the three input files
#   U1SS_RESULTS_DIR       Output directory
#   U1SS_PLUS_BIGWIG       Plus-strand U1SS bigWig file
#   U1SS_MINUS_BIGWIG      Minus-strand U1SS bigWig file
#   U1SS_GENE_TABLE        Tab-delimited gene-coordinate table
#   U1SS_NEAREST_OUTPUT    Output RData file
#   U1SS_TOP_N             Number of nearest U1 sites (default: 10)

required_packages <- c("GenomicRanges", "rtracklayer", "S4Vectors")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    "Install the following required Bioconductor package(s): ",
    paste(missing_packages, collapse = ", ")
  )
}

parse_positive_integer <- function(value, variable_name) {
  parsed_value <- suppressWarnings(as.integer(value))
  if (length(parsed_value) != 1L || is.na(parsed_value) || parsed_value < 1L) {
    stop(variable_name, " must be a positive integer; received: ", value)
  }
  parsed_value
}

data_dir <- Sys.getenv("U1SS_DATA_DIR", unset = "data")
results_dir <- Sys.getenv("U1SS_RESULTS_DIR", unset = "results")
top_n <- parse_positive_integer(
  Sys.getenv("U1SS_TOP_N", unset = "10"),
  "U1SS_TOP_N"
)

plus_bigwig <- Sys.getenv(
  "U1SS_PLUS_BIGWIG",
  unset = file.path(data_dir, "map_5ss_hg38_single_nt_plus_corrected.bw")
)
minus_bigwig <- Sys.getenv(
  "U1SS_MINUS_BIGWIG",
  unset = file.path(data_dir, "map_5ss_hg38_single_nt_minus_corrected.bw")
)
gene_table_file <- Sys.getenv(
  "U1SS_GENE_TABLE",
  unset = file.path(data_dir, "all_protein_coding_genes_MAIN_Select.tsv")
)
output_file <- Sys.getenv(
  "U1SS_NEAREST_OUTPUT",
  unset = file.path(
    results_dir,
    sprintf("nearest_%d_U1_sites_by_gene.RData", top_n)
  )
)

input_files <- c(
  plus_bigwig = plus_bigwig,
  minus_bigwig = minus_bigwig,
  gene_table = gene_table_file
)
missing_files <- input_files[!file.exists(input_files)]
if (length(missing_files) > 0L) {
  stop(
    "The following input file(s) were not found: ",
    paste(missing_files, collapse = ", ")
  )
}

dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

message("Importing strand-specific U1SS score tracks.")
u1ss_plus <- rtracklayer::import(plus_bigwig, format = "BigWig")
u1ss_minus <- rtracklayer::import(minus_bigwig, format = "BigWig")

GenomicRanges::strand(u1ss_plus) <- "+"
GenomicRanges::strand(u1ss_minus) <- "-"
u1ss_all_sites <- c(u1ss_plus, u1ss_minus)

if (length(u1ss_all_sites) == 0L) {
  stop("The U1SS bigWig files contain no intervals.")
}
if (!"score" %in% names(S4Vectors::mcols(u1ss_all_sites))) {
  stop("The imported U1SS intervals must contain a numeric `score` column.")
}

u1_scores <- suppressWarnings(
  as.numeric(S4Vectors::mcols(u1ss_all_sites)$score)
)
if (any(!is.finite(u1_scores))) {
  stop("The U1SS score tracks contain missing or non-numeric scores.")
}
S4Vectors::mcols(u1ss_all_sites)$score <- u1_scores

message("Reading protein-coding gene coordinates.")
gene_table <- utils::read.delim(
  gene_table_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required_gene_columns <- c("chr", "TSS", "strand", "gene_name")
if (!all(required_gene_columns %in% names(gene_table))) {
  stop(
    "The gene-coordinate table must contain the columns: ",
    paste(required_gene_columns, collapse = ", ")
  )
}

gene_table <- gene_table[
  !is.na(gene_table$gene_name) & nzchar(gene_table$gene_name),
  ,
  drop = FALSE
]

if (nrow(gene_table) == 0L) {
  stop("No genes with non-empty gene symbols were found.")
}
if (anyDuplicated(gene_table$gene_name)) {
  duplicated_names <- unique(gene_table$gene_name[
    duplicated(gene_table$gene_name) |
      duplicated(gene_table$gene_name, fromLast = TRUE)
  ])
  stop(
    "Gene symbols must be unique. Duplicated symbol(s): ",
    paste(utils::head(duplicated_names, 10L), collapse = ", ")
  )
}
if (any(!gene_table$strand %in% c("+", "-"))) {
  stop("The gene-coordinate table contains strand values other than `+` or `-`.")
}

gene_table$TSS <- suppressWarnings(as.integer(gene_table$TSS))
if (any(is.na(gene_table$TSS)) || any(gene_table$TSS < 1L)) {
  stop("All TSS coordinates must be positive integers in 1-based coordinates.")
}

tss_ranges <- GenomicRanges::makeGRangesFromDataFrame(
  gene_table,
  seqnames.field = "chr",
  start.field = "TSS",
  end.field = "TSS",
  strand.field = "strand"
)
names(tss_ranges) <- gene_table$gene_name

gene_loci <- paste0(
  as.character(GenomicRanges::seqnames(tss_ranges)),
  ":",
  as.character(GenomicRanges::strand(tss_ranges))
)
site_loci <- paste0(
  as.character(GenomicRanges::seqnames(u1ss_all_sites)),
  ":",
  as.character(GenomicRanges::strand(u1ss_all_sites))
)
unmatched_loci <- setdiff(unique(gene_loci), unique(site_loci))
if (length(unmatched_loci) > 0L) {
  stop(
    "No strand-matched U1SS intervals were found for these chromosome/strand ",
    "combinations: ",
    paste(unmatched_loci, collapse = ", ")
  )
}

message(
  "Finding the ", top_n,
  " nearest strand-matched U1SS intervals for ",
  length(tss_ranges), " genes."
)

nearest_indices <- GenomicRanges::nearestKNeighbors(
  tss_ranges,
  u1ss_all_sites,
  k = top_n,
  ignore.strand = FALSE
)

nearest_U1sitelist <- lapply(seq_along(tss_ranges), function(i) {
  indices <- as.integer(nearest_indices[[i]])
  indices <- indices[!is.na(indices)]

  if (length(indices) < top_n) {
    stop(
      "Only ", length(indices), " strand-matched U1SS intervals were found for ",
      names(tss_ranges)[[i]], "; expected at least ", top_n, "."
    )
  }

  gene_tss <- tss_ranges[i]
  hits <- u1ss_all_sites[indices]
  distances <- GenomicRanges::distance(
    gene_tss,
    hits,
    ignore.strand = FALSE
  )

  # Resolve equal-distance ties reproducibly using the original track order.
  hit_order <- order(distances, indices)
  hits <- hits[hit_order]
  distances <- distances[hit_order]

  data.frame(
    gene = rep(names(gene_tss), top_n),
    tss_chr = rep(as.character(GenomicRanges::seqnames(gene_tss)), top_n),
    tss_pos = rep(GenomicRanges::start(gene_tss), top_n),
    tss_strand = rep(as.character(GenomicRanges::strand(gene_tss)), top_n),
    u1_chr = as.character(GenomicRanges::seqnames(hits)),
    u1_start = GenomicRanges::start(hits),
    u1_end = GenomicRanges::end(hits),
    u1_strand = as.character(GenomicRanges::strand(hits)),
    u1_score = S4Vectors::mcols(hits)$score,
    distance = as.integer(distances),
    rank = seq_len(top_n),
    stringsAsFactors = FALSE
  )
})

analysis_parameters <- list(
  top_n = top_n,
  distance_definition = "GenomicRanges::distance with strand matching",
  plus_bigwig = normalizePath(plus_bigwig),
  minus_bigwig = normalizePath(minus_bigwig),
  gene_table = normalizePath(gene_table_file),
  created_with = R.version.string
)

save(
  nearest_U1sitelist,
  analysis_parameters,
  file = output_file,
  version = 3
)

message("Saved results to: ", normalizePath(output_file))
