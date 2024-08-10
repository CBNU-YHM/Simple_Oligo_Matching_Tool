import os
import pandas as pd

# Function to find the intersection of DNA record IDs across multiple CSV files.
def find_intersection(csv_files):
    # Initialize the intersection set with record IDs from the first CSV file.
    intersection = set(pd.read_csv(csv_files[0])['DNA record_id'])
    # Iterate over the rest of the CSV files and update the intersection set.
    for csv_file in csv_files[1:]:
        data = set(pd.read_csv(csv_file)['DNA record_id'])
        intersection.intersection_update(data)
    return intersection

# Function to analyze and save the results of intersection analysis for groups of CSV files.
def process_intersection_groups(base_directory, groups, output_directory):
    print("Processing intersection groups...")
    # Create a directory to store the intersection results if it does not exist.
    intersection_dir = os.path.join(output_directory, 'intersection')
    if not os.path.exists(intersection_dir):
        os.makedirs(intersection_dir)
        print(f"Created intersection directory: {intersection_dir}")

    # Process each group of CSV files.
    for group_name, csv_files in groups.items():
        print(f"Processing group {group_name}...")
        # Convert relative file paths to absolute paths.
        csv_files = [os.path.join(base_directory, csv_file) for csv_file in csv_files]
        # Find the intersection of DNA record IDs in the group.
        intersection = find_intersection(csv_files)
        print(f"Intersection for group {group_name}: {len(intersection)} records found.")

        # Save the intersection results to a CSV file.
        result_df = pd.DataFrame({'DNA record_id': list(intersection)})
        result_file = os.path.join(
            intersection_dir,
            f'{group_name}_intersection_{len(intersection)}_records.csv'  # Filename includes group name and record count.
        )
        result_df.to_csv(result_file, index=False)
        print(f"Intersection results saved to {result_file}.")

# Example usage of processing the intersection of DNA records from CSV files, categorized by group.
output_directory = r'C:\Users\yhm\Desktop\Data Result\Inclusivity\재계산\14. BA.2.86(Omicron) 68 sequences.fasta_Filtered_by_5%(60_68)'

# Define groups of CSV files associated with specific testing protocols or institutions.
groups = {
    'China CDC_ORF1ab': [
        'ORF1ab-F_matching.csv',
        'ORF1ab-R_matching.csv',
        'ORF1ab-P_matching.csv',
    ],
    # Additional groups can be defined similarly.
}

# Call the function to process each group and find intersections.
process_intersection_groups(output_directory, groups, output_directory)
