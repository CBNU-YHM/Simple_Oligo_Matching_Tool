import os
import csv
import concurrent.futures
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import nt_search
from functools import partial

# Function to find the best match for a primer within a given sequence.
def find_best_match(seq, primer):
    best_match = None
    best_similarity = 0
    best_position = None
    # Iterate through the sequence to find the best match based on similarity score.
    for i in range(len(seq) - len(primer) + 1):
        match_seq = seq[i:i+len(primer)]
        similarity = calculate_similarity(primer, match_seq)
        if similarity > best_similarity:
            best_match = match_seq
            best_similarity = similarity
            best_position = i
    return best_match, best_position, best_similarity

# Function to calculate the similarity percentage using IUPAC codes.
def calculate_similarity(primer, sequence):
    iupac_code = {
        'A': re.compile('A'), 'C': re.compile('C'), 'G': re.compile('G'), 'T': re.compile('T'),
        'R': re.compile('[AG]'), 'Y': re.compile('[CT]'), 'S': re.compile('[GC]'), 'W': re.compile('[AT]'),
        'K': re.compile('[GT]'), 'M': re.compile('[AC]'), 'B': re.compile('[CGT]'), 'D': re.compile('[AGT]'),
        'H': re.compile('[ACT]'), 'V': re.compile('[ACG]'), 'N': re.compile('[ACGT]')
    }
    matches = sum(bool(iupac_code[p].match(s)) for p, s in zip(primer, sequence))
    similarity = (matches / len(primer)) * 100
    return round(similarity, 2)

# Function to write headers for CSV files for matching and non-matching results.
def write_csv_headers(matching_writer, nonmatching_writer):
    matching_writer.writerow(['DNA record_id', 'primer_Sequence', 'Match_position sequence', 'Similarity'])
    nonmatching_writer.writerow(['DNA record_id', 'primer_Sequence', 'Similar_DNA_match_position', 'Similar_DNA_Sequence', 'Similarity'])
    print("CSV headers written.")

# Function to process and save results where the primer matches the DNA sequence.
def process_matching_results(record, matches, primer_sequence, matching_writer):
    similarity = "100%"
    for position in matches[1:]:
        matching_writer.writerow([record.id, primer_sequence, position + 1, similarity])

# Function to handle results where the primer does not fully match the sequence.
def process_nonmatching_results(record, primer_sequence, nonmatching_writer, match_seq, position, similarity):
    if match_seq:
        nonmatching_writer.writerow([record.id, primer_sequence, position + 1, match_seq, f"{similarity:.2f}%"])
    else:
        nonmatching_writer.writerow([record.id, primer_sequence, 'N/A', 'N/A', "0%"])

# Main processing function for each input file, writing results to CSV files.
def process_file(input_file, output_directory, primers):
    filename = os.path.basename(input_file)
    filename_without_extension = os.path.splitext(filename)[0]
    file_output_dir = os.path.join(output_directory, filename_without_extension)
    if not os.path.exists(file_output_dir):
        os.makedirs(file_output_dir)
    print(f"Processing file: {filename}")
    print(f"Created/Found directory: {file_output_dir}")

    sequence_ids = []
    records = list(SeqIO.parse(input_file, "fasta"))  # Load all records into memory beforehand.
    for record in records:
        sequence_ids.append(record.id)

    for primer_name, primer_sequence in primers.items():
        matching_csv_path = os.path.join(file_output_dir, f"{primer_name}_matching.csv")
        nonmatching_csv_path = os.path.join(file_output_dir, f"{primer_name}_nonmatching.csv")

        with open(matching_csv_path, 'w', newline='') as matching_file, \
             open(nonmatching_csv_path, 'w', newline='') as nonmatching_file:
            matching_writer = csv.writer(matching_file)
            nonmatching_writer = csv.writer(nonmatching_file)
            write_csv_headers(matching_writer, nonmatching_writer)

            for record in records:
                matches = nt_search(str(record.seq), str(primer_sequence))
                if len(matches) > 1:
                    process_matching_results(record, matches, str(primer_sequence), matching_writer)
                else:
                    match_seq, position, similarity = find_best_match(str(record.seq), str(primer_sequence))
                    if similarity == 100:
                        matching_writer.writerow([record.id, primer_sequence, position + 1, f"{similarity:.2f}%"])
                    else:
                        process_nonmatching_results(record, str(primer_sequence), nonmatching_writer, match_seq, position, similarity)

# Function to run the analysis in parallel across multiple FASTA files.
def analyze_fasta_with_primers_parallel(input_directory, output_directory, primers):
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
    print(f"Using output directory: {output_directory}")

    input_files = [os.path.join(input_directory, filename) for filename in os.listdir(input_directory) if filename.endswith(".FASTA") or filename.endswith(".fasta")]
    print(f"Found input files: {input_files}")
    with concurrent.futures.ProcessPoolExecutor(max_workers=os.cpu_count()) as executor:  # Create workers equal to the number of CPU cores.
        executor.map(partial(process_file, output_directory=output_directory, primers=primers), input_files)

if __name__ == '__main__':
    input_directory = r"C:\Users\yhm\Desktop\test"
    output_directory = r"C:\Users\yhm\Desktop\test\result0621"
    primers = {
        # Define primers with their corresponding sequences.
        "ORF1ab-F": Seq("CCCTGTGGGTTTTACACTTAA"),
        # Additional primers with their sequences go here.
    }
    analyze_fasta_with_primers_parallel(input_directory, output_directory, primers)
