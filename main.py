


import Bio
from Bio import Entrez
from Bio import Seq
from Bio import SeqIO

import ftplib

import gzip

import os
import os.path

from shutil import copyfileobj


# TODO: Make it so that it doesn't download all the data files all the time.

def main():
    get_data()

    protein_count = 0

    filenames = os.listdir("dna-sequences")
    for fn in filenames:
        # Ignore files that aren't genbank files. (For example, the compressed
        # versions of those files)
        if os.path.splitext(fn)[1] != ".gb":
            continue

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

    retrieve_files(files, "dna-sequences")

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
            with gzip.open(os.path.join(dest_dir, filename), "rb") as compressed_file:
                with open(os.path.join(dest_dir, path), "wb") as decompressed_file:
                    copyfileobj(compressed_file, decompressed_file)
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
