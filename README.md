#This tutorial was followed: https://github.com/IndexThePlanet/Logan/blob/main/Chickens.md

#getting-a-list-of-accessions - this was from the nci website (https://www.ncbi.nlm.nih.gov/sra/)  by searching candida albicans [Organsim]

#then selecting "Source --> DNA sequences" and clicking "Send results to Run Selector":

#Downloaded the following accession list associated with DNA sequences. The associated list of accessions were then used to download contigs from the LOGAN database via:

process_accessions_8.py accession_list.txt reference_ncbi.fa .

#The resulting dna sequences were then convereted to amino acid .fasta files via:
python dna_2_protein_recursive2.py .

Resulting datasets availabe: https://doi.org/10.5281/zenodo.15142125
