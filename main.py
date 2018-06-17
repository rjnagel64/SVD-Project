


import Bio
from Bio import Seq
from Bio import SeqIO
import os


def main():
    protein_count = 0

    filenames = os.listdir("dna-sequences")
    for fn in filenames:
        for record in translate_sequences("dna-sequences/" + fn):
            for protein in extract_proteins(record):
                # If there is an 'X' in a protein sequence, that means the DNA sequence
                # had some 'N's in it. Because of ambiguity, we discard the rest of the
                # proteins from this record.
                if "X" in protein:
                    break

                # TODO: Do something with the protein sequences.
                # For now, I am just counting how many there are.
                protein_count += 1

    print(f"protein_count: {protein_count} proteins")


def get_data():
    pass

def translate_sequences(filename):
    """Translate the DNA sequences in a file.

    arguments:
    - filename: The file containing the DNA sequences
    
    returns:
    An iterator for the translated sequences.
    """
    for (i, record) in enumerate(SeqIO.parse(filename, "fasta")):
        yield record.translate()

def extract_proteins(record):
    """Extract the protein sequences from a translated record.
    
    arguments:
    - record: A translated amino acid sequence

    returns:
    An iterator for each protein sequence in the record.
    """
    index = 0
    while index != -1 and index < len(record.seq):
        start_index = record.seq.find("M", index)
        if start_index == -1:
            # If we didn't find a start codon, there are no more start codons,
            # and we are done.
            break

        stop_index = record.seq.find("*", start_index)
        if stop_index == -1:
            seq = record.seq[start_index:]
        else:
            # Plus one so that the protein sequence includes the stop codon.
            seq = record.seq[start_index:stop_index + 1]

        yield seq

        index = stop_index + 1

def create_matrix():
    pass


if __name__ == "__main__":
    main()
