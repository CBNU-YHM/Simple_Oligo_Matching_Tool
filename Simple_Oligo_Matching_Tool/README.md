### Simple_Oligo_Matching_Tool V6_EN.py

- **Function**: This script analyzes how specific oligonucleotides (short DNA or RNA molecules) match within a given DNA sequence. Essentially, it compares the sequence of a given primer or oligo with the target DNA sequence to identify matching sections and provides these results to the user.
- **Usage**: Researchers can use this tool to assess the specificity of primers for certain DNA sequences or to select suitable oligos for experiments targeting specific regions.

### intersection_analysis_part V2_EN.py

- **Function**: This script, part of the Simple_Oligo_Matching_Tool, performs the task of finding the intersection of DNA record IDs from multiple CSV files. It is used to identify DNA sequences that are commonly observed across different experimental conditions or experimental groups.
- **Usage**: This analysis helps determine which DNA sequences repeatedly appear across multiple datasets, aiding in understanding the importance or ubiquity of these sequences.

### primer_verification_part V7_EN.py

- **Function**: This script, part of the Simple_Oligo_Matching_Tool, verifies how well given primers function against DNA sequences. It evaluates whether the primers accurately recognize and amplify the expected DNA regions, automating the process of experimental verification to confirm the efficacy of the primers.
- **Usage**: The script provides crucial information for optimizing or modifying primers to be used in PCR experiments, by assessing the efficiency and specificity of the primers.

### filteration_part_In_folder V1_EN.py

- **Function**: This script automates the process of filtering DNA sequences within all FASTA files in a specified directory. It evaluates sequences based on their content of a specific symbol, determined by whether they are protein or nucleotide sequences, and categorizes them into filtered and not filtered groups based on a given threshold.
- **Usage**: This tool is essential for researchers dealing with bulk sequence data, facilitating rapid segregation of sequences that meet specific criteria across entire folders. It is particularly useful in projects where preliminary data cleaning and organization are needed before more detailed analysis.

### filteration_part_Single V1_EN.py

- **Function**: This script focuses on filtering sequences from a single FASTA file, applying specific criteria to identify and categorize sequences that contain a certain percentage of a specified symbol ('X' for unknown amino acids in proteins or 'N' for any nucleotide in DNA). Sequences are then sorted into filtered and not filtered based on a predefined threshold.
- **Usage**: Ideal for detailed examination of individual sequence files, this script supports researchers in performing targeted filtration tasks. It is particularly valuable in studies requiring precise sequence editing and preparation before further genetic analysis or experimental procedures.

