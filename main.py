


import Bio
from Bio import Seq
from Bio import SeqIO
import os


def main():
    print("Hello, World!")

    translate_sequences()


def get_data():
    pass

def translate_sequences():
    filenames = os.listdir("dna-sequences")
    for fn in filenames:
        for (i, record) in enumerate(SeqIO.parse(f"dna-sequences/{fn}", "fasta")):
            SeqIO.write(record.translate(), f"aa-sequences/{i}-{fn}", "fasta")

def extract_proteins():
    filenames = os.listdir("aa-sequences")
    for fn in filenames:
        record = next(SeqIO.parse(f"aa-sequences/{fn}", "fasta"))

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

            # TODO: Record this sequence somewhere.

            index = stop_index + 1

def create_matrix():
    pass


if __name__ == "__main__":
    main()
