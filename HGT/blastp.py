import os
import sys
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO

# Check if the correct number of arguments are provided
if len(sys.argv) != 2:
    print("Usage: python script_name.py <input_fasta_file>")
    sys.exit(1)

# Get the input FASTA file from command line arguments
input_fasta_file = sys.argv[1]

# Read the FASTA file
with open(input_fasta_file) as fasta_handle:
    fasta_sequence = fasta_handle.read()

# Run the BLASTP search online against the nr database with hitlist_size set to 250
result_handle = NCBIWWW.qblast("blastp", "nr", fasta_sequence, hitlist_size=250)

# Get the base name of the input file and use it to name the output file
base_name = os.path.splitext(os.path.basename(input_fasta_file))[0]
output_file_name = base_name + ".txt"

# Create the output directory if it doesn't exist
output_dir = os.path.join(os.path.dirname(input_fasta_file), "blastp_files")
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Create the full path for the output file
output_file_path = os.path.join(output_dir, output_file_name)

# Parse the BLASTP results and write to a text file directly
with open(output_file_path, "w") as out_file:
    out_file.write("# BLASTP 2.10.1+\n")
    out_file.write("# Fields: query acc., subject acc., evalue, bit score, alignment length, % identity\n")
    out_file.write("# 250 hits found\n")

    blast_records = NCBIXML.parse(result_handle)

    try:
        for blast_record in blast_records:
            if not hasattr(blast_record, 'alignments'):
                continue
            query_acc = blast_record.query.split()[0]
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    # Calculate percent identity
                    percent_identity = (hsp.identities / hsp.align_length) * 100
                    out_file.write(f"{query_acc}\t{alignment.accession}\t{hsp.expect}\t{hsp.bits}\t{hsp.align_length}\t{percent_identity:.2f}\n")
    except ValueError as e:
        print(f"Error parsing BLAST record: {e}")

    out_file.write("# BLAST processed 1 queries\n")

result_handle.close()
