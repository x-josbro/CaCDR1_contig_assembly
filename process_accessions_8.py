import os
import sys
import subprocess
import logging
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Set up logging
logging.basicConfig(filename='processing_issues.log', level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')

def convert_sam_to_bam(sam_file, output_dir):
    accession = os.path.basename(sam_file).split('.')[0]
    bam_file = os.path.join(output_dir, f"{accession}.bam")
    sorted_bam_file = os.path.join(output_dir, f"{accession}.sorted.bam")
    
    logging.info(f"Converting SAM to BAM, sorting and indexing for {accession}...")
    
    try:
        # Convert SAM to BAM
        subprocess.run(f"samtools view -b {sam_file} > {bam_file}", shell=True, check=True)
        
        # Sort BAM
        subprocess.run(f"samtools sort {bam_file} -o {sorted_bam_file}", shell=True, check=True)
        
        # Index sorted BAM
        subprocess.run(f"samtools index {sorted_bam_file}", shell=True, check=True)
        
        # Remove the intermediate BAM file
        os.remove(bam_file)
        
        logging.info(f"Finished converting, sorting, and indexing for {accession}")
        return sorted_bam_file
    except subprocess.CalledProcessError as e:
        logging.error(f"Error in SAM to BAM conversion for {accession}: {str(e)}")
        return None

def run_samtools_consensus(bam_file, output_file):
    try:
        command = f"samtools consensus -a -aa --show-ins no --min-BQ 30 --min-MQ 40 -d 3 {bam_file} -o {output_file}"
        subprocess.run(command, shell=True, check=True)
        return True
    except subprocess.CalledProcessError as e:
        logging.error(f"Error in consensus generation for {os.path.basename(bam_file)}: {str(e)}")
        return False

def pad_consensus_sequence(reference_file, consensus_file):
    try:
        reference_record = next(SeqIO.parse(reference_file, "fasta"))
        consensus_record = next(SeqIO.parse(consensus_file, "fasta"))

        reference_seq = str(reference_record.seq)
        consensus_seq = str(consensus_record.seq)

        if len(consensus_seq) < len(reference_seq):
            padding_length = len(reference_seq) - len(consensus_seq)
            consensus_seq += 'N' * padding_length
            logging.info(f"Padded consensus sequence with {padding_length} 'N's to match reference length.")
        
        # Save the padded consensus sequence back to file
        padded_consensus_record = SeqRecord(
            Seq(consensus_seq),
            id=consensus_record.id,
            name=consensus_record.name,
            description="Padded consensus sequence"
        )
        SeqIO.write(padded_consensus_record, consensus_file, "fasta")
        return True
    except Exception as e:
        logging.error(f"Error in padding consensus sequence: {str(e)}")
        return False

def concatenate_sequences(reference_file, consensus_file, input_sequences_file):
    try:
        with open(input_sequences_file, 'w') as outfile:
            for fname in [reference_file, consensus_file]:
                with open(fname) as infile:
                    outfile.write(infile.read())
        logging.info(f"Concatenated {reference_file} and {consensus_file} into {input_sequences_file}")
        return True
    except Exception as e:
        logging.error(f"Error concatenating sequences: {str(e)}")
        return False

def align_sequences_strictly(input_sequences_file, aligned_file):
    try:
        # Align the sequences using MUSCLE
        muscle_command = [
            "muscle",
            "-in", input_sequences_file,
            "-out", aligned_file,
            "-maxiters", "32", # Higher iterations for stricter alignment
            "-diags" # Fast diagonal approximation
        ]
        subprocess.run(muscle_command, check=True)
        logging.info(f"Alignment completed strictly. Output saved to {aligned_file}")
        return True
    except subprocess.CalledProcessError as e:
        logging.error(f"Error in strict sequence alignment using MUSCLE: {str(e)}")
        return False

def fill_gaps_with_reference(aligned_file, output_file):
    try:
        # Parse the aligned sequences
        sequences = list(SeqIO.parse(aligned_file, "fasta"))
        
        if len(sequences) < 2:
            logging.error(f"Expected two sequences in the alignment file: {aligned_file}")
            return False
        
        reference_seq = str(sequences[0].seq)
        consensus_seq = str(sequences[1].seq)

        logging.info(f"Before gap filling - Consensus length: {len(consensus_seq)}, Reference length: {len(reference_seq)}")
        logging.info(f"Number of 'N's in consensus before filling: {consensus_seq.count('N')}")

        # Initialize the list for the filled sequence
        filled_seq = []

        # Fill gaps/Ns in the consensus sequence with reference sequence, avoiding unexpected insertions
        for i, (cons, ref) in enumerate(zip(consensus_seq, reference_seq)):
            if cons == 'N' or cons == '-':
                filled_seq.append(ref)
            else:
                filled_seq.append(cons)
        
        filled_seq = ''.join(filled_seq)

        # Check if the filled sequence is the same length as the reference and consensus
        if len(filled_seq) != len(reference_seq):
            logging.error(f"Filled sequence length {len(filled_seq)} does not match reference length {len(reference_seq)}")

        logging.info(f"After gap filling - Filled sequence length: {len(filled_seq)}")
        logging.info(f"Number of filled positions: {sum(1 for cons, ref in zip(consensus_seq, filled_seq) if cons == 'N' or cons == '-')}")

        filled_record = SeqRecord(
            Seq(filled_seq),
            id=sequences[1].id,
            name=sequences[1].name,
            description="Consensus sequence with gaps filled from reference"
        )
        SeqIO.write(filled_record, output_file, "fasta")
        return True
    except Exception as e:
        logging.error(f"Error in gap filling for {aligned_file}: {str(e)}")
        return False

def process_accession(accession, reference_file, root_dir):
    for dirpath, dirnames, filenames in os.walk(root_dir):
        sam_file = os.path.join(dirpath, f"{accession}.minimap2_output")
        if os.path.exists(sam_file):
            logging.info(f"Processing {accession}...")
            
            # Convert SAM to sorted BAM
            sorted_bam_file = convert_sam_to_bam(sam_file, dirpath)
            if not sorted_bam_file:
                logging.error(f"Failed to process BAM file for {accession}")
                return False

            consensus_file = os.path.join(dirpath, f"{accession}_consensus.fa")
            filled_consensus_file = os.path.join(dirpath, f"{accession}_consensus_filled.fa")
            aligned_file = os.path.join(dirpath, f"{accession}_aligned_sequences.fas")
            input_sequences_file = os.path.join(dirpath, f"{accession}_input_sequences.fasta")

            if not run_samtools_consensus(sorted_bam_file, consensus_file):
                logging.error(f"Failed to generate consensus for {accession}")
                return False

            # Pad the consensus sequence to match the reference length
            if not pad_consensus_sequence(reference_file, consensus_file):
                logging.error(f"Failed to pad consensus sequence for {accession}")
                return False

            # Concatenate reference and consensus sequences
            if not concatenate_sequences(reference_file, consensus_file, input_sequences_file):
                logging.error(f"Failed to concatenate sequences for {accession}")
                return False

            # Align sequences strictly using MUSCLE
            if not align_sequences_strictly(input_sequences_file, aligned_file):
                logging.error(f"Failed to strictly align sequences for {accession}")
                return False

            # Correct 'N's and gaps in the consensus sequence using the reference
            if not fill_gaps_with_reference(aligned_file, filled_consensus_file):
                logging.error(f"Failed to fill gaps for {accession}")
                return False

            logging.info(f"Finished processing {accession}")
            return True
    
    logging.warning(f"Could not find .minimap2_output file for {accession}")
    return False

def main(accession_list_file, reference_file, root_dir):
    with open(accession_list_file, 'r') as f:
        accessions = [line.strip() for line in f]
    
    total = len(accessions)
    successful = 0
    failed = 0

    for i, accession in enumerate(accessions, 1):
        logging.info(f"Processing accession {i} of {total}: {accession}")
        if process_accession(accession, reference_file, root_dir):
            successful += 1
        else:
            failed += 1
    
    logging.info(f"Processing complete. Successful: {successful}, Failed: {failed}, Total: {total}")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py <accession_list_file> <reference_file> <root_directory>")
        sys.exit(1)
    
    accession_list_file = sys.argv[1]
    reference_file = sys.argv[2]
    root_dir = sys.argv[3]
    
    main(accession_list_file, reference_file, root_dir)
