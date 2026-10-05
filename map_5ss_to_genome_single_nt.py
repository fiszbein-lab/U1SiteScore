import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from concurrent.futures import ProcessPoolExecutor, as_completed
import logging

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(message)s')
logger = logging.getLogger()

# Load the reference genome
def load_genome(genome_fasta):
    genome = {}
    for record in SeqIO.parse(genome_fasta, "fasta"):
        if "GL" not in record.id and "KI" not in record.id:  # Filter out chromosomes with "GL" or "KI"
            genome[record.id] = str(record.seq)
    return genome

# Function to calculate the window around a given point
def get_window_sequence(genome, chrom, position, strand='+'):
    if strand == '+':
        start = max(0, position - 4)  # Shift -1
        end = position + 5
    else:  # strand == '-'
        start = max(0, position - 0)  # Shift +3
        end = position + 9

    seq = genome[chrom][start:end]
    if strand == '-':
        seq = str(Seq(seq).reverse_complement())
    return seq
    
# Function to map 5' splice site scores and calculate sum within the window
def map_splice_scores(genome, chrom, position, scores_dict, strand='+'):
    kmer = get_window_sequence(genome, chrom, position, strand)
    score = scores_dict.get(kmer, 0)  # Get score or 0 if kmer not in dict
    return score

# Load 9-mer sequences and their scores
def load_scores(scores_file):
    scores_dict = {}
    with open(scores_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            kmer = parts[0]
            score = float(parts[1])
            scores_dict[kmer] = score
    return scores_dict

# Function to process a whole chromosome and condense the BEDGRAPH data
def process_chromosome(genome, scores_dict, chrom, step_size=1, strand='+'):
    chrom_length = len(genome[chrom])
    segment_data = []
    current_start = 0
    current_score = map_splice_scores(genome, chrom, current_start, scores_dict, strand)
    current_end = current_start + step_size

    for position in range(step_size, chrom_length - 6, step_size):
        total_score = map_splice_scores(genome, chrom, position, scores_dict, strand)

        if total_score == current_score:
            current_end = position + step_size
        else:
            segment_data.append((chrom, current_start, current_end, current_score))
            current_start = position
            current_end = position + step_size
            current_score = total_score

    # Append the last segment if not empty
    segment_data.append((chrom, current_start, current_end, current_score))

    # Condense segments further if needed
    condensed_data = []
    for i, (chrom, start, end, score) in enumerate(segment_data):
        if i == 0:
            condensed_data.append((chrom, start, end, score))
        else:
            last_chrom, last_start, last_end, last_score = condensed_data[-1]
            if chrom == last_chrom and score == last_score and start == last_end:
                # Extend the previous segment
                condensed_data[-1] = (last_chrom, last_start, end, score)
            else:
                condensed_data.append((chrom, start, end, score))

    return condensed_data

# Write results to BEDGRAPH format
def write_bedgraph(output_file, bedgraph_data):
    with open(output_file, 'w') as f:
        f.write("track type=bedGraph\n")
        for chrom, start, end, score in bedgraph_data:
            f.write(f"{chrom}\t{start}\t{end}\t{score:.2f}\n")  # Format score to two decimal places

# Main function to process the whole genome
def process_genome(genome_fasta, scores_file, output_prefix, step_size=1, num_workers=4):
    genome = load_genome(genome_fasta)
    scores_dict = load_scores(scores_file)
    plus_bedgraph_data = []
    minus_bedgraph_data = []

    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        futures = []
        for chrom in genome:
            logger.info(f"Processing chromosome: {chrom}")
            futures.append(executor.submit(process_chromosome, genome, scores_dict, chrom, step_size, '+'))
            futures.append(executor.submit(process_chromosome, genome, scores_dict, chrom, step_size, '-'))

        for future in as_completed(futures):
            result = future.result()
            if futures.index(future) % 2 == 0:
                plus_bedgraph_data.extend(result)
            else:
                minus_bedgraph_data.extend(result)

    # Ensure data is sorted by chromosome and position
    plus_bedgraph_data.sort(key=lambda x: (x[0], x[1]))
    minus_bedgraph_data.sort(key=lambda x: (x[0], x[1]))

    write_bedgraph(f"{output_prefix}_plus.bg", plus_bedgraph_data)
    write_bedgraph(f"{output_prefix}_minus.bg", minus_bedgraph_data)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Map 5' splice site scores to the genome and export in BEDGRAPH format.")
    parser.add_argument("--fasta", required=True, help="Path to the reference genome FASTA file.")
    parser.add_argument("--scores", required=True, help="Path to the file with 9-mer sequences and their scores.")
    parser.add_argument("--outname", required=True, help="Output file prefix for BEDGRAPH files.")
    parser.add_argument("--num_threads", type=int, default=4, help="Number of worker processes to use.")
    parser.add_argument("--bin", type=int, default=1, help="Step size for sliding window (should be 1 for single nucleotide resolution).")
    
    args = parser.parse_args()

    process_genome(args.fasta, args.scores, args.outname, args.bin, args.num_threads)
