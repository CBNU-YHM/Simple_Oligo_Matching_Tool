# Simple_Oligo_Matching_Tool
Simple_Oligo_Matching_Tool

This manual provides a guide to using the provided Python script for analyzing DNA sequences with specific primers. The script leverages the BioPython library for sequence processing, concurrent programming for handling multiple files in parallel, and pandas for data manipulation.

### Required Library Installations

Before running the script, you need to install the necessary Python libraries. You can install these using pip, Python's package installer. Open your command prompt or terminal and enter the following commands:

- **BioPython** (for DNA sequence processing):
  ```bash
  pip install biopython
  ```
- **Pandas** (for data handling):
  ```bash
  pip install pandas
  ```

### File Storage Locations

When you run the analysis script, it processes DNA sequences from the files you specify and generates output in the form of CSV files. These files are automatically stored in a structured way so you can easily find and manage the results.

- **Output Directory**: Before you start, you'll specify an output directory. This is a folder on your computer where all the results will be saved.
- **Folder Structure**: For each FASTA file processed, the script creates a new subfolder in the output directory. The subfolder is named after the FASTA file, minus its file extension. This helps keep results organized, especially when processing multiple files.
- **Files in Each Folder**: Inside each subfolder, you'll find:
  - A CSV file for sequences that perfectly match the primer.
  - A CSV file for sequences that don't perfectly match but have the best partial match.
  - A CSV file listing all the sequence IDs analyzed from the FASTA file.

**Example**: If your output directory is set to `C:\DNA\Results` and you're analyzing `sample.fasta`, you will find the results in `C:\DNA\Results\sample\`.

### Primer Input

Primers are short DNA sequences used to locate a specific sequence within a DNA sample. In this script, you'll need to specify which primers you want to use for the analysis.

- **Defining Primers**: You define primers in the script as a dictionary, where each key-value pair consists of a primer name and its sequence. This allows the script to know which sequences to look for in your DNA samples.
- **How to Input Primers**: In the script, there's a section where you can list all the primers you want to use. You simply add entries to the dictionary in the format `"PrimerName": "PrimerSequence"`.

**Example**:
```python
primers = {
    "Primer1": "ATCG",
    "Primer2": "GCTA"
}
```
This tells the script that `Primer1` has the sequence `ATCG`, and `Primer2` has the sequence `GCTA`.

### User Defined Paths and Settings

This part of the script allows you to customize which directories (folders) and specific settings to use for the analysis. It's how you tell the script where to find the DNA files and how to organize the output.

- **Defining Paths and Settings**: You provide these details in a list of dictionaries. Each dictionary specifies a base path for the DNA files, along with a set of primers relevant to that path.
- **Structure of the Dictionary**:
  - `base_path`: This is the folder where your FASTA files are located.
  - `primer_sets`: A dictionary that maps a descriptive name to a list of primer names you want to use for files in this base path.

**Example**:
```python
user_defined_paths_settings = [
    {
        "base_path": "path/to/your/FASTA/files",
        "primer_sets": {
            'Set1': ['Primer1', 'Primer2'],
        }
    },
    # Additional configurations can be added here
]
```
In this example, `base_path` is where the script looks for FASTA files, and `primer_sets` is which primers to apply to these files. `Set1` is just a label; it could be anything descriptive.

### Tips for Ease of Use

- **Prepare Your Files**: Ensure all FASTA files are in the specified `base_path` before running the script.
- **Check Primer Sequences**: Make sure the primer sequences you enter are correct to ensure accurate results.
- **Run from the Correct Location**: Open your command prompt or terminal in the script's directory, or navigate to it using `cd path/to/your/script`.

By setting up everything as described, you can smoothly run the script and get your results organized and accessible, making it easier to conduct your analyses.


### Detailed Function Descriptions

- **find_best_match(seq, primer)**:
  - **Purpose**: Identifies the best local match for a primer within a DNA sequence.
  - **Parameters**:
    - `seq`: The DNA sequence to search within.
    - `primer`: The primer sequence to match against the DNA sequence.
  - **Returns**: A tuple containing the best matching sequence, its position, and the similarity percentage.

- **calculate_similarity(seq1, seq2)**:
  - **Purpose**: Calculates the similarity percentage between two sequences based on IUPAC codes.
  - **Parameters**:
    - `seq1`, `seq2`: The sequences to compare.
  - **Returns**: The similarity percentage as a float.

- **write_csv_headers(matching_writer, nonmatching_writer)**:
  - **Purpose**: Writes headers to the CSV files for matching and non-matching results.
  - **Parameters**:
    - `matching_writer`: CSV writer for the matching results file.
    - `nonmatching_writer`: CSV writer for the non-matching results file.

- **process_matching_results(record, matches, primer_sequence, matching_writer)**:
  - **Purpose**: Processes and writes results for DNA sequences that fully match the primer.
  - **Parameters**:
    - `record`: The DNA record being analyzed.
    - `matches`: Positions where the primer fully matches.
    - `primer_sequence`: The primer sequence used.
    - `matching_writer`: CSV writer for the matching results.

- **process_nonmatching_results(record, primer_sequence, nonmatching_writer, match_seq, position, similarity)**:
  - **Purpose**: Handles the logging of non-full matches.
  - **Parameters**:
    - `record`: The DNA record being analyzed.
    - `primer_sequence`: The primer used.
    - `nonmatching_writer`: CSV writer for non-matching results.
    - `match_seq`: The best matching sequence.
    - `position`: Position of the best match.
    - `similarity`: Similarity of the best match.

- **analyze_fasta_with_primers_parallel(input_directory, output_directory, primers)**:
  - **Purpose**: Analyzes all FASTA files in a directory with specified primers.
  - **Parameters**:
    - `input_directory`: Path to the directory containing FASTA files.
    - `output_directory`: Path where output results will be stored.
    - `primers`: Dictionary of primers to use for analysis.

- **find_intersection(csv_files)**:
  - **Purpose**: Finds the intersection of DNA record IDs across multiple CSV files.
  - **Parameters**:
    - `csv_files`: List of CSV files to analyze.
  - **Returns**: A set of common DNA record IDs.

- **process_intersection_groups(base_directory, groups, output_directory)**:
  - **Purpose**: Processes groups of CSV files to find common DNA records.
  - **Parameters**:
    - `base_directory`: Base directory for the analysis.
    - `groups`: Dictionary of groups and their associated CSV files.
    - `output_directory`: Directory to store intersection results.

  * groups: Dictionary of groups and their associated CSV files.
  * output_directory: Directory to store intersection results.



  
