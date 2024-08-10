from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import nt_search
import os
import csv
import re
import pandas as pd
from pathlib import Path
import concurrent.futures
from functools import partial

# Function to find the best local match for a given primer within a DNA sequence.
# It iterates through the sequence to find the highest similarity match for the primer.
def find_best_match(seq, primer):
    best_match = None
    best_similarity = 0
    best_position = None
    # Loop through the sequence to find the best match for the primer.
    for i in range(len(seq) - len(primer) + 1):
        match_seq = seq[i:i+len(primer)]
        similarity = calculate_similarity(primer, match_seq)
        if similarity > best_similarity:
            best_match = match_seq
            best_similarity = similarity
            best_position = i
    # Return the best matching sequence, its position, and similarity percentage.
    if best_match:
        return best_match, best_position, best_similarity
    else:
        return None, None, 0

# Calculates the similarity percentage between two sequences based on IUPAC nucleotide codes.
def calculate_similarity(seq1, seq2):
    iupac_code = {
        'A': re.compile('A'), 'C': re.compile('C'), 'G': re.compile('G'), 'T': re.compile('T'),
        'R': re.compile('[AG]'), 'Y': re.compile('[CT]'), 'S': re.compile('[GC]'), 'W': re.compile('[AT]'),
        'K': re.compile('[GT]'), 'M': re.compile('[AC]'), 'B': re.compile('[CGT]'), 'D': re.compile('[AGT]'),
        'H': re.compile('[ACT]'), 'V': re.compile('[ACG]'), 'N': re.compile('[ACGT]')
    }
    matches = sum(bool(iupac_code[p].match(s)) for p, s in zip(seq1, seq2))
    return (matches / len(seq1)) * 100

# Write headers for both matching and non-matching results to CSV files.
def write_csv_headers(matching_writer, nonmatching_writer):
    matching_writer.writerow(['DNA record_id', 'primer_Sequence', 'Match_position sequence', 'Similarity'])
    nonmatching_writer.writerow(['DNA record_id', 'primer_Sequence', 'Similar_DNA_match_position', 'Similar_DNA_Sequence', 'Similarity'])

# Process and record results for DNA sequences that fully match the primer.
def process_matching_results(record, matches, primer_sequence, matching_writer):
    similarity = "100%"
    for position in matches[1:]:
        matching_writer.writerow([record.id, primer_sequence, position + 1, similarity])

# Record results for DNA sequences that do not completely match the primer.
# Includes best partial matches with their respective positions and similarity scores.
def process_nonmatching_results(record, primer_sequence, nonmatching_writer, match_seq, position, similarity):
    if match_seq:
        nonmatching_writer.writerow([record.id, primer_sequence, position + 1, match_seq, f"{similarity:.2f}%"])
    else:
        nonmatching_writer.writerow([record.id, primer_sequence, 'N/A', 'N/A', "0%"])

# Analyzes all FASTA files in a specified directory using the given primers, and writes the results to CSV files.
def analyze_fasta_with_primers_parallel(input_directory, output_directory, primers):
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
        print(f"Created output directory: {output_directory}")

    input_files = [os.path.join(input_directory, filename) for filename in os.listdir(input_directory)
                   if filename.endswith(".FASTA") or filename.endswith(".fasta")]

    # Utilizes a process pool to handle multiple files concurrently.
    with concurrent.futures.ProcessPoolExecutor() as executor:
        executor.map(partial(process_file, output_directory=output_directory, primers=primers), input_files)

# Finds the intersection of DNA record IDs across multiple CSV files to identify common sequences.
def find_intersection(csv_files):
    intersection = set(pd.read_csv(csv_files[0])['DNA record_id'])
    for csv_file in csv_files[1:]:
        data = set(pd.read_csv(csv_file)['DNA record_id'])
        intersection.intersection_update(data)
    return intersection

# Analyzes and records the results of intersection analysis for groups of CSV files.
def process_intersection_groups(base_directory, groups, output_directory):
    print("Processing intersection groups...")
    for group_name, csv_files in groups.items():
        print(f"Processing group {group_name}...")
        intersection = find_intersection(csv_files)
        print(f"Intersection for group {group_name}: {len(intersection)} records found.")

        result_df = pd.DataFrame({'DNA record_id': list(intersection)})
        intersection_dir = Path(output_directory) / 'intersection' # / group_name
        intersection_dir.mkdir(parents=True, exist_ok=True)

        result_file = intersection_dir / f'{group_name}_intersection_{len(result_df)}count.csv'
        result_df.to_csv(result_file, index=False)
        print(f"Intersection results for {group_name} saved to {result_file}.")

# Generates a dictionary of groups dynamically based on user-defined path settings for output directory management.
def build_groups_dict(user_defined_paths_settings, output_directory):
    groups = {}
    for setting in user_defined_paths_settings:
        base_path = setting["base_path"]
        for set_name, primers in setting["primer_sets"].items():
            group_key = f"{set_name}_{base_path.replace(' ', '_').replace('%', 'Perc')}"
            groups[group_key] = [
                str(Path(output_directory) / base_path / f"{primer}_matching.csv")
                for primer in primers
            ]
    return groups

def process_file(input_file, output_directory, primers):
    filename = os.path.basename(input_file)
    print(f"Processing file: {filename}")
    filename_without_extension = os.path.splitext(filename)[0]
    file_output_dir = os.path.join(output_directory, filename_without_extension)
    if not os.path.exists(file_output_dir):
        os.makedirs(file_output_dir)
        print(f"Created directory for {filename}: {file_output_dir}")

    sequence_ids = []
    for record in SeqIO.parse(input_file, "fasta"):
        sequence_ids.append(record.id)

    for primer_name, primer_sequence in primers.items():
        print(f"Analyzing {filename} with primer {primer_name}...")
        matching_csv_path = os.path.join(file_output_dir, f"{primer_name}_matching.csv")
        nonmatching_csv_path = os.path.join(file_output_dir, f"{primer_name}_nonmatching.csv")

        with open(matching_csv_path, 'w', newline='') as matching_file, \
             open(nonmatching_csv_path, 'w', newline='') as nonmatching_file:
            matching_writer = csv.writer(matching_file)
            nonmatching_writer = csv.writer(nonmatching_file)
            write_csv_headers(matching_writer, nonmatching_writer)

            for record in SeqIO.parse(input_file, "fasta"):
                matches = nt_search(str(record.seq), str(primer_sequence))
                if len(matches) > 1:
                    process_matching_results(record, matches, str(primer_sequence), matching_writer)
                else:
                    match_seq, position, similarity = find_best_match(str(record.seq), str(primer_sequence))
                    process_nonmatching_results(record, str(primer_sequence), nonmatching_writer, match_seq, position, similarity)

    # Create a CSV file listing all sequence IDs processed.
    sequence_id_file = os.path.join(file_output_dir, f"{filename_without_extension}_{len(sequence_ids)}count.csv")
    with open(sequence_id_file, 'w', newline='') as seq_id_file:
        writer = csv.writer(seq_id_file)
        writer.writerow(['Sequence ID'])
        for seq_id in sequence_ids:
            writer.writerow([seq_id])
    print(f"Created a CSV file with sequence IDs: {sequence_id_file}")

    print(f"Finished analyzing {filename} with primer {primer_name}.")

# Example usage: configure paths and primers for running the analysis.
if __name__ == '__main__':
    user_defined_paths_settings = [
        {
            "base_path": "2024-03-13 Beta B.1.351 1189 sequences.fasta_Filtered_by_5%(1118_1189)",
            "primer_sets": {
                'NIID_2019-nCOV_N': ['NIID_2019-nCOV_N_F2', 'NIID_2019-nCOV_N_R2', 'NIID_2019-nCOV_N_P2'],
            }
        },
        {
            "base_path": "2024-03-13 Beta B.1.351 1189 sequences.fasta_Filtered_by_5%(1118_1189) - 복사본",
            "primer_sets": {
                'NIID_2019-nCOV_N': ['NIID_2019-nCOV_N_F2', 'NIID_2019-nCOV_N_R2', 'NIID_2019-nCOV_N_P2'],
            }
        },
        # Additional path and primer set configurations can be added here.
    ]

    input_directory = r"C:\Users\yhm\Desktop\test\2024-08-04 test"
    output_directory = r"C:\Users\yhm\Desktop\test\2024-08-04 test"

    # Generate a dictionary of groups based on user-defined path settings.
    groups = build_groups_dict(user_defined_paths_settings, output_directory)

    # Set up primers.
    primers = {
        "NIID_2019-nCOV_N_F2": Seq("AAATTTTGGGGACCAGGAAC"),
        "NIID_2019-nCOV_N_R2": Seq("GTTGACCTACACAGCTGCCA"),
        "NIID_2019-nCOV_N_P2": Seq("ATGTCGCGCATTGGCATGGA"),
    }

    # Analyze FASTA files and process the results.
    analyze_fasta_with_primers_parallel(input_directory, output_directory, primers)
    process_intersection_groups(output_directory, groups, output_directory)
