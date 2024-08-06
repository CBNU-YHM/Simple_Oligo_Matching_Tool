# Simple_Oligo_Matching_Tool

### Introduction to the User Manual

This manual is designed to facilitate users in utilizing a specialized Python script for analyzing DNA sequences with specific primers. The script makes extensive use of the BioPython library for efficient sequence processing, employs concurrent programming techniques to manage multiple files simultaneously, and utilizes pandas for sophisticated data manipulation.

The manual is divided into two main sections:

1. **Script(Simple_Oligo_Matching_Tool) Usage Guide**: Provides detailed instructions on how to use the script, including descriptions of its functions, setup of input parameters, and how to interpret the outputs.
2. **IDE Setup Instructions**: Offers a comprehensive guide to installing and configuring Visual Studio Code (VS Code), a popular and accessible Integrated Development Environment (IDE), to run the Python script effectively. This section is tailored to help even those with no previous experience in programming or development environments.

Each section is structured to provide step-by-step instructions to ensure ease of use and to help users efficiently achieve their objectives with the script.

## Script(Simple_Oligo_Matching_Tool) Usage Guide

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



## Introduction to Setting Up an Integrated Development Environment(IDE) Setup Instructions

This section of the manual is dedicated to assisting users in establishing a functional Integrated Development Environment (IDE) to effectively run Python scripts. We have chosen Visual Studio Code (VS Code) as the recommended IDE due to its broad accessibility and user-friendly interface. This guide aims to provide step-by-step instructions to ensure even users with no prior experience in programming environments can set up and begin using VS Code for Python development with ease.

### Step 1: Install Python

Before using any Python scripts, you need to have Python installed on your computer:

1. **Download Python**:
   - Visit the official Python website at [python.org](https://www.python.org/downloads/).
   - Download the latest version of Python for your operating system (Windows, macOS, or Linux).
   - Follow the installation instructions on the website.

2. **Check Installation**:
   - Open your command prompt (Windows) or terminal (macOS or Linux).
   - Type `python --version` and press Enter. If Python is installed correctly, you should see the version number displayed.

### Step 2: Install Visual Studio Code (VS Code)

VS Code is a popular IDE that makes coding easier with its intuitive interface and powerful features:

1. **Download VS Code**:
   - Go to the [VS Code website](https://code.visualstudio.com/) and download the installer for your operating system.
   - Run the installer and follow the on-screen instructions to complete the installation.

2. **Open VS Code**:
   - Once installed, open VS Code from your applications or programs list.

### Step 3: Set Up VS Code for Python Development

To make VS Code ready for Python development, you need to install the Python extension:

1. **Install Python Extension**:
   - In VS Code, go to the Extensions view by clicking on the square icon on the sidebar or pressing `Ctrl+Shift+X`.
   - Search for "Python" and find the extension provided by Microsoft.
   - Click on "Install" to add the extension to VS Code.

2. **Open Your Project**:
   - In VS Code, go to `File > Open Folder...` and select the folder where you've saved your Python script.

### Step 4: Configure the Python Interpreter

Ensure that VS Code uses the correct Python interpreter:

1. **Open the Command Palette**:
   - Press `Ctrl+Shift+P` to open the Command Palette.
   - Type "Python: Select Interpreter" and press Enter.

2. **Select Python Interpreter**:
   - Choose the Python version you installed earlier. It usually appears with the path where Python is installed.

### Step 5: Run the Script

Now that everything is set up, you can run the Python script:

1. **Open the Script**:
   - In VS Code, navigate to the Explorer on the sidebar, and double-click your Python script (it should end in `.py`) to open it.

2. **Run the Script**:
   - Right-click in the editor window containing your script and select "Run Python File in Terminal" or use the play button in the top-right corner of the editor.
   - VS Code will execute the Python script in the terminal at the bottom of the window, and you'll see the outputs or any messages printed by your script.

### Tips for Success

- **Check Your Script's Requirements**: Make sure you install any required libraries (like `biopython` or `pandas`) using pip in the VS Code terminal. For example, you can type `pip install biopython` directly in the terminal.
- **Save Your Changes**: Remember to save your script after any modifications by pressing `Ctrl+S`.

By following these steps, even users with no prior knowledge of development environments can set up and run a Python script efficiently, allowing them to focus more on their research or tasks.
  
