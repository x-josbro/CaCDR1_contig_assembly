import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Seq import translate

def find_longest_orf(sequence):
    longest_orf = ""
    for strand, nuc in [(+1, sequence), (-1, sequence.reverse_complement())]:
        for frame in range(3):
            for pro in nuc[frame:].translate(table=1).split("*"):
                if len(pro) > len(longest_orf):
                    longest_orf = pro
    return longest_orf

def process_fasta(input_file, output_file):
    with open(input_file, "r") as input_handle, open(output_file, "w") as output_handle:
        for record in SeqIO.parse(input_handle, "fasta"):
            sequence = record.seq

            # Trim the sequence to a multiple of three if necessary
            if len(sequence) % 3 != 0:
                trim_len = len(sequence) % 3
                sequence = sequence[:-trim_len]
                print(f"Trimmed {trim_len} bases from the end of {record.id} to make sequence length a multiple of three.")

            longest_orf = find_longest_orf(sequence)
            protein_record = SeqRecord(
                Seq(longest_orf),
                id=record.id + "_protein",
                description="Longest ORF from " + record.description
            )
            SeqIO.write(protein_record, output_handle, "fasta")

def convert_dna_to_protein_recursive(root_dir):
    for dirpath, dirnames, filenames in os.walk(root_dir):
        for filename in filenames:
            if filename.endswith("_consensus_filled.fa"):
                input_file = os.path.join(dirpath, filename)
                output_file = os.path.join(dirpath, filename.replace(".fa", "_aa.fa"))
                print(f"Converting {input_file} to {output_file}")
                try:
                    process_fasta(input_file, output_file)
                    print(f"Successfully converted {input_file}")
                except Exception as e:
                    print(f"Error processing {input_file}: {str(e)}")

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 2:
        print("Usage: python script.py <root_directory>")
        sys.exit(1)
    
    root_directory = sys.argv[1]
    convert_dna_to_protein_recursive(root_directory)
    print("Conversion process completed.")
