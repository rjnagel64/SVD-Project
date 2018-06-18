


import Bio
from Bio import Seq
from Bio import SeqIO

import ftplib

import gzip

import os
import os.path

from shutil import copyfileobj

# TODO: One of the sequences in the drosophila melanogaster genome (sequence 0)
# has a length that is not a multiple of three. Biopython suggests appending an
# 'N' to the end of the sequence to rectify this.

# 1869 sequences for drosophila melanogaster.
# 63 sequences for thalassiosira pseudonana. Hmm. (But significantly larger on average)

DNA_SEQUENCE_DIR = "dna-sequences"
PROTEIN_SEQUENCE_DIR = "protein-sequences"

def main(args):
    # Only download the files if we don't already have them.
    if not os.path.exists(DNA_SEQUENCE_DIR):
        get_data()

    write_proteins()


def get_data():
    """Retrieve data from the NCBI servers and databanks."""

    files = [
        # Drosophila melanogaster is a fruit fly.
        ("genomes/all/GCF/000/001/215/GCF_000001215.4_Release_6_plus_ISO1_MT/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.gbff.gz",
        "drosophila_melanogaster.gb.gz"),
        # Thalassiosira pseudonana is some sort of diatom, with a genome that's
        # somewhat smaller than Drosophila, but it has more chromosomes.
        ("genomes/all/GCF/000/149/405/GCF_000149405.2_ASM14940v2/GCF_000149405.2_ASM14940v2_genomic.gbff.gz",
        "thalassiosira_pseudonana.gb.gz"),
    ]

    retrieve_files(files, DNA_SEQUENCE_DIR)

# TODO: Error handling if a file does not exist on the server?
def retrieve_files(filepaths, dest_dir):
    """Retrieve a number of files from GenBank via FTP.
    
    arguments:
    - filepaths: A list of `(file_path, local_name)` pairs, where `file_path` is
        the path to the file on the server, and `local_name` is the name that the
        file will have once it is downloaded.
    - dest_dir: A local path to the place where the files will be downloaded.
        It will be created if it does not already exist.
    """

    # Create the destination directory, if necessary.
    if not os.path.exists(dest_dir):
        os.mkdir(dest_dir)

    f = ftplib.FTP("ftp.ncbi.nlm.nih.gov")
    f.login()

    for (fp, fn) in filepaths:
        (path, filename) = os.path.split(fp)
        print(f"Downloading {filename} ...")
        # Return to the root directory of the available files.
        f.cwd("/")
        # Go to the directory containing the file we want.
        f.cwd(path)
        with open(os.path.join(dest_dir, fn), "wb") as fd:
            f.retrbinary("RETR " + filename, fd.write)
    
    f.quit()

    print("Files downloaded.")

    # Decompress the files downloaded.
    for filename in os.listdir(dest_dir):
        (path, ext) = os.path.splitext(filename)
        if ext == ".gz":
            print(f"Decompressing {filename}...")
            compressed_filename = os.path.join(dest_dir, filename)
            
            # Decompress the file and copy its contents into a new file.
            with gzip.open(compressed_filename, "rb") as compressed_file:
                with open(os.path.join(dest_dir, path), "wb") as decompressed_file:
                    copyfileobj(compressed_file, decompressed_file)

            # Delete the compressed file. It's not useful.
            os.remove(compressed_filename)

        else:
            # Ignore files that are not compressed.
            continue

    print("Files downloaded and decompressed.")

def translate_sequences(filename):
    """Translate the DNA sequences in a file.

    arguments:
    - filename: The file containing the DNA sequences, in GenBank format
    
    returns:
    An iterator for the translated sequences.
    """
    for (i, record) in enumerate(SeqIO.parse(filename, "genbank")):
        print(f"Translating sequence {i} of {filename}...")
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
            # If a sequence ends without the stop codon, it's not a protein.
            return
        else:
            # Plus one so that the protein sequence includes the stop codon.
            seq = str(record.seq[start_index:stop_index + 1])

        yield seq

        index = stop_index + 1

def write_proteins():
    """Write the proteins from a given file into a file."""
    if not os.path.exists(PROTEIN_SEQUENCE_DIR):
        os.mkdir(PROTEIN_SEQUENCE_DIR)

    total_proteins = 0
    num_sequences = 0

    for filename in os.listdir(DNA_SEQUENCE_DIR):
        (basename, ext) = os.path.splitext(filename)

        # Ignore files not in the genbank format.
        if ext != ".gb":
            continue

        print(f"Processing file {filename}...")

        with open(os.path.join(PROTEIN_SEQUENCE_DIR, basename), "w") as f:
            for (i, record) in enumerate(translate_sequences(os.path.join(DNA_SEQUENCE_DIR, filename))):
                print(f"Extracting proteins from sequence number {i}...")
                num_proteins = 0
                num_sequences += 1

                for protein in extract_proteins(record):
                    num_proteins += 1
                    f.write(protein + "\n")

                print(f"\t{num_proteins} proteins in that sequence.")
                total_proteins += num_proteins

    print(f"There were {total_proteins} proteins in total, over {num_sequences} sequences.")

def create_matrix():
    pass


if __name__ == "__main__":
    from sys import argv
    main(argv)
