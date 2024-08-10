from Bio import SeqIO
import os

def create_directories(base_path, threshold):
    # Create directories for storing filtered and not filtered sequences based on the threshold percentage
    # Paths are constructed using absolute paths
    filtered_path = os.path.join(base_path, f'result/filtered/filtered_{threshold}%')
    not_filtered_path = os.path.join(base_path, f'result/filtered/not_filtered_{threshold}%')
    os.makedirs(filtered_path, exist_ok=True)  # Ensure directories are created if they do not exist
    os.makedirs(not_filtered_path, exist_ok=True)
    return filtered_path, not_filtered_path

def get_symbol_for_sequence_type(sequence_type):
    # Determine the symbol to count based on the type of sequence: protein or nucleotide
    if sequence_type.lower() == 'protein':
        return 'X'  # 'X' used for unknown amino acids in protein sequences
    elif sequence_type.lower() == 'nucleotide':
        return 'N'  # 'N' used for any nucleotide in DNA sequences
    else:
        raise ValueError("Invalid sequence type. Please enter 'protein' or 'nucleotide'.")

def filter_sequences(input_file, sequence_type, threshold):
    base_path = os.path.dirname(input_file)
    filtered_path, not_filtered_path = create_directories(base_path, threshold)

    symbol = get_symbol_for_sequence_type(sequence_type)  # Retrieve the appropriate symbol for counting

    # Temporary file paths for filtered and not filtered sequences
    filtered_file_temp_path = os.path.join(filtered_path, f"Filtered_temp.fasta")
    not_filtered_file_temp_path = os.path.join(not_filtered_path, f"Not_Filtered_temp.fasta")

    # Initialize counters for filtered, not filtered, and total sequences processed
    filtered_count = 0
    not_filtered_count = 0
    total_count = 0

    # Process sequences and apply filtering based on the threshold
    with open(filtered_file_temp_path, 'w') as filtered_file, open(not_filtered_file_temp_path, 'w') as not_filtered_file:
        for record in SeqIO.parse(input_file, 'fasta'):
            sequence = str(record.seq).upper()
            n_count = sequence.count(symbol)  # Count the occurrences of the symbol
            n_ratio = (n_count / len(sequence)) * 100  # Calculate the percentage of the symbol
            total_count += 1

            if n_ratio < threshold:  # If the symbol ratio is below the threshold, consider it filtered
                SeqIO.write(record, filtered_file, 'fasta')
                filtered_count += 1
            else:  # If the symbol ratio meets or exceeds the threshold, it is not filtered
                SeqIO.write(record, not_filtered_file, 'fasta')
                not_filtered_count += 1

    # Rename temporary files to their final names with additional information
    filtered_file_final_path = os.path.join(filtered_path, f"{os.path.basename(input_file).replace('.fasta', '')}_Filtered_by_{threshold}%({filtered_count}_{total_count}).fasta")
    not_filtered_file_final_path = os.path.join(not_filtered_path, f"{os.path.basename(input_file).replace('.fasta', '')}_Not_Filtered_by_{threshold}%({not_filtered_count}_{total_count}).fasta")

    os.rename(filtered_file_temp_path, filtered_file_final_path)
    os.rename(not_filtered_file_temp_path, not_filtered_file_final_path)

    print(f"Filtered sequences: {filtered_file_final_path}")
    print(f"Not filtered sequences: {not_filtered_file_final_path}")

# User configuration settings
input_path = r"C:\Users\ADMIN\Desktop\SARS-CoV-2 Seq(2024-03)\2024-03-13 omicron parent B.1.1.529 2056 sequences.fasta"
sequence_type = "Nucleotide"  # User specifies either 'Nucleotide' or 'Protein'
n_threshold = 5  # Threshold percentage for filtering based on 'N' count

# Function call
filter_sequences(input_path, sequence_type, n_threshold)
