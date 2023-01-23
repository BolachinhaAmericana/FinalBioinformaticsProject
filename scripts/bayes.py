import os

# Set the input and output directories
input_dir = "/path/to/input/fasta/files"
output_dir = "/path/to/output/tree/files"

# Loop through all fasta files in the input directory
for file in os.listdir(input_dir):
    if file.endswith(".fasta"):
        # Get the file name without the extension
        file_name = os.path.splitext(file)[0]
        # Run Mr bayes on the fasta file
        os.system(f"mrbayes -s {input_dir}/{file} -n {file_name} -o {output_dir}/{file_name}.tre")

# Print a message to indicate the script has completed
print("Mr bayes tree creation complete!")