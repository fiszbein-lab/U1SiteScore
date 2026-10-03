#!/usr/bin/env python3
"""Map 9-mer 5' splice-site scores across a reference genome.

For positions sampled at a fixed step size, the script sums scores for all
9-mers in a strand-specific window and writes separate plus- and minus-strand
bedGraph files. Adjacent intervals with identical scores are merged.

Example
-------
python map_5ss_to_genome_publication.py \
    --fasta data/hg38.fa \
    --scores data/u1ss_9mer_scores.tsv \
    --output-prefix results/map_5ss_hg38 \
    --workers 8
"""

import argparse
from collections import OrderedDict
from concurrent.futures import ProcessPoolExecutor, as_completed
import logging
import math
import os
from typing import Dict, Iterable, List, Sequence, Tuple

try:
    from Bio import SeqIO
    from Bio.Seq import Seq
except ImportError as exc:
    raise SystemExit(
        "Biopython is required. Install it with `python -m pip install biopython`."
    ) from exc


KMER_LENGTH = 9
DEFAULT_EXCLUDED_CONTIG_MARKERS = ("GL", "KI")
BedGraphRecord = Tuple[str, int, int, float]

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
)
LOGGER = logging.getLogger(__name__)


def load_genome(
    genome_fasta: str,
    include_alternate_contigs: bool = False,
) -> "OrderedDict[str, str]":
    """Read FASTA sequences while preserving their input order.

    By default, contig identifiers containing ``GL`` or ``KI`` are excluded to
    reproduce the primary-contig filtering in the original analysis.
    """
    genome = OrderedDict()

    with open(genome_fasta, "r") as fasta_handle:
        for record in SeqIO.parse(fasta_handle, "fasta"):
            if (
                not include_alternate_contigs
                and any(marker in record.id for marker in DEFAULT_EXCLUDED_CONTIG_MARKERS)
            ):
                continue
            if record.id in genome:
                raise ValueError(
                    "Duplicate FASTA sequence identifier: {0}".format(record.id)
                )
            genome[record.id] = str(record.seq).upper()

    if not genome:
        raise ValueError("No reference sequences remained after FASTA filtering.")

    return genome


def load_scores(scores_file: str) -> Dict[str, float]:
    """Read a two-column, tab-delimited table of 9-mers and numeric scores."""
    scores = {}

    with open(scores_file, "r") as score_handle:
        for line_number, raw_line in enumerate(score_handle, start=1):
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue

            fields = line.split("\t")
            if len(fields) != 2:
                raise ValueError(
                    "Expected two tab-delimited fields at line {0} of {1}.".format(
                        line_number, scores_file
                    )
                )

            kmer = fields[0].upper()
            if len(kmer) != KMER_LENGTH:
                raise ValueError(
                    "Expected a {0}-mer at line {1}; found {2!r}.".format(
                        KMER_LENGTH, line_number, fields[0]
                    )
                )
            if kmer in scores:
                raise ValueError(
                    "Duplicate k-mer {0!r} at line {1}.".format(kmer, line_number)
                )

            try:
                score = float(fields[1])
            except ValueError as exc:
                raise ValueError(
                    "Non-numeric score at line {0}: {1!r}.".format(
                        line_number, fields[1]
                    )
                ) from exc

            if not math.isfinite(score):
                raise ValueError(
                    "Non-finite score at line {0}: {1!r}.".format(
                        line_number, fields[1]
                    )
                )
            scores[kmer] = score

    if not scores:
        raise ValueError("The score table contains no 9-mer scores.")

    return scores


def get_window_sequence(
    chromosome_sequence: str,
    position: int,
    window_size: int,
    strand: str,
) -> str:
    """Return a bounded, strand-oriented sequence around a zero-based position."""
    start = max(0, position - window_size)
    end = min(len(chromosome_sequence), position + window_size)
    sequence = chromosome_sequence[start:end]

    if strand == "-":
        return str(Seq(sequence).reverse_complement())
    if strand != "+":
        raise ValueError("strand must be `+` or `-`")
    return sequence


def map_splice_scores(
    chromosome_sequence: str,
    position: int,
    scores: Dict[str, float],
    window_size: int,
    strand: str,
) -> float:
    """Sum 9-mer scores in one strand-oriented sequence window."""
    window_sequence = get_window_sequence(
        chromosome_sequence,
        position,
        window_size,
        strand,
    )
    return sum(
        scores.get(window_sequence[index : index + KMER_LENGTH], 0.0)
        for index in range(max(0, len(window_sequence) - KMER_LENGTH + 1))
    )


def process_chromosome(
    chromosome: str,
    chromosome_sequence: str,
    scores: Dict[str, float],
    window_size: int,
    step_size: int,
    strand: str,
) -> List[BedGraphRecord]:
    """Calculate and merge sampled window scores for one chromosome and strand."""
    chromosome_length = len(chromosome_sequence)
    if chromosome_length == 0:
        return []

    current_start = 0
    current_end = min(step_size, chromosome_length)
    current_score = map_splice_scores(
        chromosome_sequence,
        current_start,
        scores,
        window_size,
        strand,
    )
    records = []

    for position in range(step_size, chromosome_length, step_size):
        score = map_splice_scores(
            chromosome_sequence,
            position,
            scores,
            window_size,
            strand,
        )
        interval_end = min(position + step_size, chromosome_length)

        if score == current_score and position == current_end:
            current_end = interval_end
        else:
            records.append(
                (chromosome, current_start, current_end, current_score)
            )
            current_start = position
            current_end = interval_end
            current_score = score

    records.append((chromosome, current_start, current_end, current_score))
    return records


def process_chromosome_pair(
    chromosome: str,
    chromosome_sequence: str,
    scores: Dict[str, float],
    window_size: int,
    step_size: int,
) -> Tuple[str, List[BedGraphRecord], List[BedGraphRecord]]:
    """Process both strands while transferring a chromosome only once to a worker."""
    plus_records = process_chromosome(
        chromosome,
        chromosome_sequence,
        scores,
        window_size,
        step_size,
        "+",
    )
    minus_records = process_chromosome(
        chromosome,
        chromosome_sequence,
        scores,
        window_size,
        step_size,
        "-",
    )
    return chromosome, plus_records, minus_records


def iter_records_in_fasta_order(
    records_by_chromosome: Dict[str, List[BedGraphRecord]],
    chromosome_order: Sequence[str],
) -> Iterable[BedGraphRecord]:
    """Yield records in FASTA contig order and ascending coordinate order."""
    for chromosome in chromosome_order:
        for record in records_by_chromosome[chromosome]:
            yield record


def write_bedgraph(
    output_file: str,
    records: Iterable[BedGraphRecord],
    score_precision: int,
    include_track_line: bool,
) -> None:
    """Write validated, zero-based half-open intervals in bedGraph format."""
    with open(output_file, "w") as output_handle:
        if include_track_line:
            output_handle.write("track type=bedGraph\n")
        for chromosome, start, end, score in records:
            output_handle.write(
                "{0}\t{1}\t{2}\t{3:.{precision}f}\n".format(
                    chromosome,
                    start,
                    end,
                    score,
                    precision=score_precision,
                )
            )


def process_genome(
    genome_fasta: str,
    scores_file: str,
    output_prefix: str,
    window_size: int = 100,
    step_size: int = 10,
    workers: int = 4,
    score_precision: int = 2,
    include_alternate_contigs: bool = False,
    include_track_line: bool = True,
) -> Tuple[str, str]:
    """Run the genome-wide analysis and return the two output paths."""
    if window_size < 1:
        raise ValueError("window_size must be a positive integer")
    if step_size < 1:
        raise ValueError("step_size must be a positive integer")
    if workers < 1:
        raise ValueError("workers must be a positive integer")
    if score_precision < 0:
        raise ValueError("score_precision must be zero or greater")

    LOGGER.info("Loading reference genome: %s", genome_fasta)
    genome = load_genome(genome_fasta, include_alternate_contigs)
    LOGGER.info("Loaded %d reference sequences", len(genome))

    LOGGER.info("Loading 9-mer scores: %s", scores_file)
    scores = load_scores(scores_file)
    LOGGER.info("Loaded %d unique 9-mer scores", len(scores))

    chromosome_order = list(genome.keys())
    plus_results = {}
    minus_results = {}

    with ProcessPoolExecutor(max_workers=workers) as executor:
        futures = {
            executor.submit(
                process_chromosome_pair,
                chromosome,
                sequence,
                scores,
                window_size,
                step_size,
            ): chromosome
            for chromosome, sequence in genome.items()
        }

        for future in as_completed(futures):
            expected_chromosome = futures[future]
            chromosome, plus_records, minus_records = future.result()
            if chromosome != expected_chromosome:
                raise RuntimeError(
                    "Worker returned results for an unexpected chromosome."
                )
            plus_results[chromosome] = plus_records
            minus_results[chromosome] = minus_records
            LOGGER.info("Finished chromosome: %s", chromosome)

    output_directory = os.path.dirname(os.path.abspath(output_prefix))
    os.makedirs(output_directory, exist_ok=True)
    plus_output = "{0}_plus.bg".format(output_prefix)
    minus_output = "{0}_minus.bg".format(output_prefix)

    write_bedgraph(
        plus_output,
        iter_records_in_fasta_order(plus_results, chromosome_order),
        score_precision,
        include_track_line,
    )
    write_bedgraph(
        minus_output,
        iter_records_in_fasta_order(minus_results, chromosome_order),
        score_precision,
        include_track_line,
    )

    LOGGER.info("Wrote plus-strand bedGraph: %s", plus_output)
    LOGGER.info("Wrote minus-strand bedGraph: %s", minus_output)
    return plus_output, minus_output


def positive_integer(value: str) -> int:
    """Argparse converter for integers greater than zero."""
    integer = int(value)
    if integer < 1:
        raise argparse.ArgumentTypeError("value must be a positive integer")
    return integer


def nonnegative_integer(value: str) -> int:
    """Argparse converter for integers greater than or equal to zero."""
    integer = int(value)
    if integer < 0:
        raise argparse.ArgumentTypeError("value must be zero or greater")
    return integer


def build_argument_parser() -> argparse.ArgumentParser:
    """Construct the command-line interface."""
    parser = argparse.ArgumentParser(
        description=(
            "Map 9-mer 5' splice-site scores to strand-specific genomic "
            "windows and export bedGraph files."
        )
    )
    parser.add_argument(
        "--fasta",
        required=True,
        help="Reference-genome FASTA file.",
    )
    parser.add_argument(
        "--scores",
        required=True,
        help="Two-column TSV containing 9-mer sequences and scores.",
    )
    parser.add_argument(
        "--output-prefix",
        "--outname",
        dest="output_prefix",
        required=True,
        help="Output prefix; `_plus.bg` and `_minus.bg` are appended.",
    )
    parser.add_argument(
        "--workers",
        "--num-threads",
        "--num_threads",
        dest="workers",
        type=positive_integer,
        default=4,
        help="Worker processes (default: 4).",
    )
    parser.add_argument(
        "--window-size",
        "--window",
        dest="window_size",
        type=positive_integer,
        default=100,
        help="Bases on each side of a sampled position (default: 100).",
    )
    parser.add_argument(
        "--step-size",
        "--bin",
        dest="step_size",
        type=positive_integer,
        default=10,
        help="Distance in bases between sampled positions (default: 10).",
    )
    parser.add_argument(
        "--score-precision",
        type=nonnegative_integer,
        default=2,
        help="Digits after the decimal point in bedGraph scores (default: 2).",
    )
    parser.add_argument(
        "--include-alternate-contigs",
        action="store_true",
        help="Include FASTA records whose identifiers contain GL or KI.",
    )
    parser.add_argument(
        "--omit-track-line",
        action="store_true",
        help="Do not write a UCSC `track type=bedGraph` header line.",
    )
    return parser


def main() -> None:
    """Parse arguments and execute the analysis."""
    args = build_argument_parser().parse_args()
    process_genome(
        genome_fasta=args.fasta,
        scores_file=args.scores,
        output_prefix=args.output_prefix,
        window_size=args.window_size,
        step_size=args.step_size,
        workers=args.workers,
        score_precision=args.score_precision,
        include_alternate_contigs=args.include_alternate_contigs,
        include_track_line=not args.omit_track_line,
    )


if __name__ == "__main__":
    main()
