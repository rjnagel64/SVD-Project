


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
        # TODO: Do some stuff with this. Actually get the protein sequences out of it.

def create_matrix():
    pass


if __name__ == "__main__":
    main()
