


import Bio
from Bio import Seq
from Bio import SeqIO

import numpy as np

from collections import Counter

from ftplib import FTP

import gzip

import os
import os.path

from shutil import copyfileobj

# TODO: One of the sequences in the drosophila melanogaster genome (sequence 0)
# has a length that is not a multiple of three. Biopython suggests appending an
# 'N' to the end of the sequence to rectify this.

DNA_SEQUENCE_DIR = "dna-sequences"
PROTEIN_SEQUENCE_DIR = "protein-sequences"
MATRIX_SAVE_DIR = "matrix"

def main(args):
    # Only download the files if we don't already have them.
    if not os.path.exists(DNA_SEQUENCE_DIR):
        get_data()

    write_proteins()

    print("Creating matrix...")
    (mat, proteins, genomes) = create_matrix2()

    print(mat.shape)

    print("Saving data...")
    save_data(mat, proteins, genomes)

    # TODO: Do SVD of the matrix

def save_data(mat, proteins, genomes):
    """Save the term-document matrix into the subdirectory `MATRIX_SAVE_DIR`."""
    if not os.path.exists(MATRIX_SAVE_DIR):
        os.mkdir("matrix")

    np.save(os.path.join(MATRIX_SAVE_DIR, "matrix.npy"), mat)

    with open(os.path.join(MATRIX_SAVE_DIR, "proteins.txt"), "w") as f:
        f.writelines(proteins)

    with open(os.path.join(MATRIX_SAVE_DIR, "genomes.txt"), "w") as f:
        f.writelines(genomes)

def get_data():
    """Retrieve data from the NCBI servers and databanks."""

    # The data files we want have information about them stored in 'genomes.txt'.
    # Perform simple parsing of that file to extract the relevant data.
    with open("genomes.txt", "r") as f:
        genomes = f.read()

    # Take only the nonempty, noncomment lines of the file.
    lines = list(filter(lambda l: l != "" and l[0] != "#", map(str.strip, genomes.splitlines())))
    # This turns ['a', 'b', 'c', 'd', ...] into [('a', 'b'), ('c', 'd'), ...].
    files = list(zip(lines[::2], lines[1::2]))

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

    f = FTP("ftp.ncbi.nlm.nih.gov")
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
    An iterator for the translated sequences, after non-coding regions are
    removed.
    """
    for (i, record) in enumerate(SeqIO.parse(filename, "genbank")):
        print(f"Translating sequence {i} of {filename}...")

        features = [f for f in record.features if f.type in ("exon", "CDS")]
        print(f"record {i}: {len(features)} features.")

        r = Seq.MutableSeq("")

        print(f"Splicing record {i}...")
        for f in features:
            r.extend(f.extract(record.seq))

        yield r.toseq().translate()

def extract_proteins(sequence):
    """Extract the protein sequences from a translated sequence.

    arguments:
    - sequence: A translated amino acid sequence

    returns:
    An iterator for each protein sequence in the translated sequence.
    """

    index = 0
    while index != -1 and index < len(sequence):
        start_index = sequence.find("M", index)
        if start_index == -1:
            # If we didn't find a start codon, there are no more start codons,
            # and we are done.
            break

        stop_index = sequence.find("*", start_index)
        if stop_index == -1:
            # If a sequence ends without the stop codon, it's not a protein.
            return
        else:
            # Plus one so that the protein sequence includes the stop codon.
            seq = str(sequence[start_index:stop_index + 1])

        yield seq

        index = stop_index + 1

def write_proteins():
    """Write the proteins from a given file into a file."""
    if not os.path.exists(PROTEIN_SEQUENCE_DIR):
        os.mkdir(PROTEIN_SEQUENCE_DIR)

    totals = Counter()
    for filename in os.listdir(DNA_SEQUENCE_DIR):
        (basename, ext) = os.path.splitext(filename)

        # Ignore files not in the genbank format.
        if ext not in (".gb", ".gbk"):
            continue

        print(f"Processing file {filename}...")

        counter = Counter()
        with open(os.path.join(PROTEIN_SEQUENCE_DIR, basename), "w") as f:
            for (i, seq) in enumerate(translate_sequences(os.path.join(DNA_SEQUENCE_DIR, filename))):
                print(f"Extracting proteins from sequence number {i}...")

                for protein in extract_proteins(seq):
                    f.write(protein + "\n")
                    counter[protein] += 1

        num_proteins = sum(counter.values())
        num_unique = len(counter.keys())
        print(f"This file contained {num_proteins} proteins, of which {num_unique} were unique.")
        totals.update(counter)

    total_proteins = sum(totals.values())
    total_unique = len(totals.keys())
    print(f"In total, there were {total_proteins} proteins, of which {total_unique} were unique.")


# This function should probably return the matrix variable.
def create_matrix():
    # TODO: The class 'Counter' from the standard library module 'collections'
    # would do nicely for counting proteins.
    # TODO: Use proper matrix types from numpy/scipy (specifically, 'ndarray')

    files = os.listdir(PROTEIN_SEQUENCE_DIR)
    allproteins = set()

    for filename in files:
       #TODO: os.path
       filepath = os.path.join(PROTEIN_SEQUENCE_DIR,filename)
       filestream = open(filepath, "r")
       for line in filestream:
           line = line.strip()
           allproteins.add(line)

       filestream.close()

    allproteins = list(allproteins)

    matrix = []
    realsequence = []

    for filename in files:
        filepath = os.path.join(PROTEIN_SEQUENCE_DIR,filename)
        sequence = open(filepath, "r")
        for line in sequence:
           line = line.strip()
           realsequence.append(line)
           vector = []
           for x in allproteins:
               cnt = realsequence.count(x)
               vector.append(cnt)
        matrix.append(vector)
        vector = []
        realsequence = []

    print('\n'.join(['\t'.join([str(cell) for cell in row]) for row in matrix]))

def create_matrix2():
    """Create the term-document matrix.

    returns:
    A tuple `(mat, proteins, genomes)` such that:
    - `mat` is a numpy ndarray.
    - `proteins` is the list of all protein sequences such that the `i`th row
        of `mat` corresponds to `proteins[i]`.
    - `genomes` is the list of genome names such that the `j`th column of `mat`
        corresponds to `genomes[j]`."""
    # `counters` associates each filename with a `Counter` of all the proteins
    # in that file.
    counters = {}
    for filename in os.listdir(PROTEIN_SEQUENCE_DIR):
        with open(os.path.join(PROTEIN_SEQUENCE_DIR, filename), "r") as f:
            counters[filename] = Counter(f)

    # Now, get a set of all the proteins (from all the files)
    all_proteins = set()
    for (filename, counter) in counters.items():
        for protein in counter.keys():
            all_proteins.add(protein)

    proteins = []
    genomes = list(counters.keys())
    mat = np.zeros((len(all_proteins), len(counters.keys())))

    # Put the genomes on the outside loop so we have a few long iterations
    # instead of many short iterations
    for (j, (genome, counter)) in enumerate(counters.items()):
        for (i, protein) in enumerate(all_proteins):
            # If a key is not found in a counter, it defaults to zero.
            mat[i, j] = counter[protein]

            # Only build the list of proteins once.
            if j == 0:
                proteins.append(protein)

    return (mat, proteins, genomes)

if __name__ == "__main__":
    from sys import argv
    main(argv)
