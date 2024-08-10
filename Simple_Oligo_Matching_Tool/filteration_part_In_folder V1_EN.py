from Bio import SeqIO
import os

def create_directories(output_directory, threshold):
    # Create directories for filtered and not filtered sequences based on the threshold
    filtered_path = os.path.join(output_directory, f'filtered_{threshold}%')
    not_filtered_path = os.path.join(output_directory, f'not_filtered_{threshold}%')
    os.makedirs(filtered_path, exist_ok=True)  # Ensure the directory is created if it doesn't exist
    os.makedirs(not_filtered_path, exist_ok=True)
    return filtered_path, not_filtered_path

def get_symbol_for_sequence_type(sequence_type):
    # Return the appropriate symbol based on the sequence type (protein or nucleotide)
    if sequence_type.lower() == 'protein':
        return 'X'  # Placeholder for unknown amino acids
    elif sequence_type.lower() == 'nucleotide':
        return 'N'  # Placeholder for any nucleotide
    else:
        raise ValueError("Invalid sequence type. Please enter 'protein' or 'nucleotide'.")

def filter_sequences(input_directory, sequence_type, threshold):
    # Main function to filter sequences based on their content of a specific symbol
    output_directory = os.path.join(input_directory, 'result/filtered')
    filtered_path, not_filtered_path = create_directories(output_directory, threshold)
    symbol = get_symbol_for_sequence_type(sequence_type)

    for filename in os.listdir(input_directory):
        if filename.endswith(".fasta") or filename.endswith(".FASTA"):
            input_file_path = os.path.join(input_directory, filename)

            filtered_records = []
            not_filtered_records = []

            for record in SeqIO.parse(input_file_path, 'fasta'):
                sequence = str(record.seq).upper()
                symbol_count = sequence.count(symbol)
                symbol_ratio = (symbol_count / len(sequence)) * 100  # Calculate the percentage of the symbol

                if symbol_ratio < threshold:
                    filtered_records.append(record)  # Add to filtered if below threshold
                else:
                    not_filtered_records.append(record)  # Add to not filtered if above or equal to threshold

            if filtered_records:
                # Create a path and save filtered records
                filtered_file_path = os.path.join(filtered_path,
                                                  f"{filename}_Filtered_by_{threshold}%({len(filtered_records)}_{len(filtered_records) + len(not_filtered_records)}).fasta")
                SeqIO.write(filtered_records, filtered_file_path, 'fasta')
                print(f"Filtered sequences: {filtered_file_path}")

            if not_filtered_records:
                # Create a path and save not filtered records
                not_filtered_file_path = os.path.join(not_filtered_path,
                                                      f"{filename}_Not_Filtered_by_{threshold}%({len(not_filtered_records)}_{len(filtered_records) + len(not_filtered_records)}).fasta")
                SeqIO.write(not_filtered_records, not_filtered_file_path, 'fasta')
                print(f"Not filtered sequences: {not_filtered_file_path}")

# User configuration settings
input_directory = r"C:\Users\ADMIN\Desktop\Data Result\VS Primer Blast"  # Folder path
sequence_type = "Nucleotide"  # User input: 'Nucleotide' or 'Protein'
n_threshold = 5  # Threshold for the percentage of 'N' or 'X'

# Call the function
filter_sequences(input_directory, sequence_type, n_threshold)
