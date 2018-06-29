
import pdb

import Bio
from Bio import Seq
from Bio import SeqIO

import matplotlib.pyplot as plt
import numpy as np
from scipy import linalg as dense_linalg
from scipy import sparse
from scipy.sparse import linalg

from collections import Counter

from ftplib import FTP

import gzip

import os
import os.path
import glob

from shutil import copyfileobj

# TODO: Add an optional normalization step to matrix creation.
#   * Create a new function, `normalize(mat, proteins, genomes)`?
#   * Do TF-IDF normalization as one option

# TODO: Create visualization of the rank 2 approximation using matplotlib.

# TODO: One of the sequences in the drosophila melanogaster genome (sequence 0)
# has a length that is not a multiple of three. Biopython suggests appending an
# 'N' to the end of the sequence to rectify this.

DNA_SEQUENCE_DIR = "dna-sequences"
PROTEIN_SEQUENCE_DIR = "protein-sequences"
MATRIX_SAVE_DIR = "matrix"
DEFAULT_GENOME_PATHS_FILE = "genomes.txt"

def main(args):
    (mat, proteins, genomes) = get_matrix(args)

    print(mat.shape)

    (W, ss, Vh) = sparse.linalg.svds(mat, k=2, which="LM")
    print(ss)
    print(Vh)

    mat = mat
    xr = (-25, 25)
    yr = (-10, 300)
    create_plot(mat, ss, Vh, xrange=xr, yrange=yr, name="test.pdf")

    mat2 = inverse_frequency(mat)
    xr = (-50, 50)
    yr = (-10, 500)
    create_plot(mat2, ss, Vh, xrange=xr, yrange=yr, name="test2.pdf")

    mat3 = norm(mat)
    xr = (-2, 2)
    yr = (-2, 2)
    create_plot(mat3, ss, Vh, xrange=xr, yrange=yr, name="test3.pdf")

def create_plot(mat, ss, Vh, xrange=None, yrange=None, name="test.pdf"):
    """Create a scatter plot of words in the new basis

    arguments:
    - mat: The term-document matrix
    - ss: Singular values (not used? may remove this parameter)
    - Vh: The matrix on the right in SVD. It forms a basis for V
    - xrange: A tuple `(xmin, xmax)` or None. Points outside this range are discarded.
    - yrange: A tuple `(ymin, ymay)` or None. Points outside this range are discarded.
    """
    if xrange is None:
        xrange = (-float("inf"), float("inf"))

    if yrange is None:
        yrange = (-float("inf"), float("inf"))

    print(f"Creating plot {name}...")
    s0, s1 = ss[0:2]

    # The rows of Vh form a basis in the new space. I am keeping them as a matrix
    # so that I can transform points into this new basis via multiplication.
    new_basis = Vh[0:2]

    # (2 * 7) @ (N * 7).T = (2 * N)

    # TODO: Get rid of outliers.
    # There is only one point with y > 1000 (~y = 4000)
    # Almost all points lie between x = -25 and x = 25
    # Also, the massive number of points is interesting for PDF renderers to deal
    # with.

    # Transform each point with matrix multiplication.
    xs, ys = new_basis @ mat.transpose()
    # Unfortunately, matplotlib insists on the inputs to the scatter plot function
    # being two separate lists of coordinates. (Not iterators, lists. *sigh*)
    xs = list(xs)
    ys = list(ys)

    new_xs = []
    new_ys = []
    for (x, y) in zip(xs, ys):
        if x < xrange[0] or x > xrange[1] or y < yrange[0] or y > yrange[1]:
            continue
        else:
            new_xs.append(x)
            new_ys.append(y)


    (fig, ax) = plt.subplots()

    # TODO: Set axis labels
    # ax.set_xlabel(...)
    # ax.set_ylabel(...)

    ax.scatter(
        new_xs,
        new_ys,
        s=2,
    )

    # TODO: Add labels to some important points?
    # Probably not. Amino acid sequences make terrible labels.

    print(f"Saving plot...")
    plt.savefig(name, format="pdf")

def get_matrix(args):
    """Obtain the term-document matrix, either by loading NCBI FTP paths from a file, all (.gbk) paths in a subdirectory, or creating
    from scratch."""
    #If arguments are passed, then use them to determine the .gbk files to use
    if len(args)>1:
        path = args[1]
        # Check if the argument given is a valid directory
        if os.path.isdir(path):
            # If so, use this directory to find the .gbk files of interest
            print("checking "+path+" for .gbk files...")
            path_list = glob.glob(os.path.join(path,"**","*.gbk"),recursive=True)
            print("found "+str(len(path_list))+" .gbk files...")
            if not os.path.exists(PROTEIN_SEQUENCE_DIR):
                print(f"Writing proteins...")
                # pdb.set_trace()
                write_proteins(path_list)
        if os.path.isfile(path):
            # If the path points to a file, use the file to search the NCBI server
            print("Searching NCBI for genome data listed in "+ path+"...")
            get_ncbi_data(path)
    else:
        # Proceed with the default pipeline via the FTP server provided by NCBI using the default list of genomes
        print("No valid local path given. Searching NCBI for genome data listed in "+ DEFAULT_GENOME_PATHS_FILE+"...")
        # If the matrix already exists, just load it.
        if os.path.exists(MATRIX_SAVE_DIR):
            print(f"Loading matrix from save...")
            mat = sparse.load_npz(os.path.join(MATRIX_SAVE_DIR, "matrix.npz"))

            print(f"Loading proteins from save...")
            with open(os.path.join(MATRIX_SAVE_DIR, "proteins.txt"), "r") as f:
                proteins = list(f)

            print(f"Loading genomes from save...")
            #Note that the filename used here is not to be confused with the default file listing FTP paths on the NCBI server
            with open(os.path.join(MATRIX_SAVE_DIR, "genomes.txt"), "r") as f:
                genomes = list(f)

            print(f"Matrix loaded.")
            return (mat, proteins, genomes)

        # Otherwise, the matrix must be created.
        print(f"Creating term-document matrix...")

        # Only download the files if we don't already have them.
        if not os.path.exists(DNA_SEQUENCE_DIR):
            print(f"Acquiring data...")
            get_ncbi_data()

        if not os.path.exists(PROTEIN_SEQUENCE_DIR):
            print(f"Writing proteins...")
            write_proteins()

    print("Constructing matrix...")
    (mat, proteins, genomes) = create_matrix2()

    print(f"Matrix created: {mat.shape}.")

    print("Saving data...")
    save_data(mat, proteins, genomes)

    return (mat, proteins, genomes)

def save_data(mat, proteins, genomes):
    """Save the term-document matrix into the subdirectory `MATRIX_SAVE_DIR`."""
    if not os.path.exists(MATRIX_SAVE_DIR):
        os.mkdir("matrix")

    #  np.save(os.path.join(MATRIX_SAVE_DIR, "matrix.npy"), mat)
    sparse.save_npz(os.path.join(MATRIX_SAVE_DIR, "matrix.npz"), mat)

    with open(os.path.join(MATRIX_SAVE_DIR, "proteins.txt"), "w") as f:
        f.writelines(proteins)

    with open(os.path.join(MATRIX_SAVE_DIR, "genomes.txt"), "w") as f:
        f.writelines(genomes)

def get_ncbi_data(path=DEFAULT_GENOME_PATHS_FILE):
    """Retrieve data from the NCBI servers and databanks."""

    # The data files we want have information about them stored in 'genomes.txt'.
    # Perform simple parsing of that file to extract the relevant data.
    with open(path, "r") as f:
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

def write_proteins(path_list=None):
    """Write the proteins from a given list of paths to .gbk files into files in the protein sequence directory."""
    if path_list == None:
        path_list = map(lambda name: os.path.join(DNA_SEQUENCE_DIR, name), os.listdir(DNA_SEQUENCE_DIR))
    if not os.path.exists(PROTEIN_SEQUENCE_DIR):
        os.mkdir(PROTEIN_SEQUENCE_DIR)

    totals = Counter()
    for path in path_list:
        (root, ext) = os.path.splitext(path)
        filename = os.path.splitext(os.path.basename(path))[0]

        # Ignore files not in the genbank format.
        if ext not in (".gb", ".gbk"):
            continue

        print(f"Processing file {path}...")

        counter = Counter()
        with open(os.path.join(PROTEIN_SEQUENCE_DIR, filename), "w") as f:
            for (i, seq) in enumerate(translate_sequences(path)):
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
    - `mat` is a scipy sparse array.
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
    genomes = list(map(lambda x: x + "\n", counters.keys()))
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

    return (sparse.csr_matrix(mat), proteins, genomes)

def inverse_frequency(mat):
    """Apply TF-IDF to a matrix.

    arguments:
    - mat: a term-document matrix with shape `(n, m)`. `m[i, j]` is how many times
    word `i` appears in document `j`.

    returns:
    Another matrix with shape `(n, m)`, where each element is the term frequency-
    inverse document frequency of that word and document."""
    # Using the raw-count definition of term frequency, tf(t, d) is just
    # mat[t, d].

    # This compares each element of `mat` to zero, and creates a matrix where
    # `nonzero[i, j] = (mat[i, j] != 0)`.
    nonzero = mat != 0
    # This counts how many `True` values were in each row of `nonzero` -- that
    # is, how many nonzero values were in that row of `mat`.
    # `counts` is an `n x 1` matrix.
    counts = nonzero.sum(axis=1)
    # To get the inverse document frequency, the number of documents is divided
    # by each count, and then the logarithm is taken.
    # `idf` is also an `n x 1` matrix.
    idf = np.log(mat.shape[1] / counts)

    # The `multiply` method works elementwise, and expands `idf` so the dimensions
    # match. The end result (since `mat` is the term frequency matrix) is a matrix
    # of TF-IDF values.
    return mat.multiply(idf)

def norm(mat):
    """Normalize each row in `mat`

    arguments:
    - mat: a sparse term-document matrix with shape `m x n`.

    returns:
    a sparse `m x n` matrix such that each entry is the corresponding entry of
    `mat` divided by that row's norm."""

    # norms : m (not m x 1, for some reason)
    norms = linalg.norm(mat, axis=1)
    # Reshape it so that it actually is m x 1.
    norms.resize(norms.shape + (1,))

    # Use elementwise division to divide each row by its norm, and remember to
    # keep the matrix sparse.
    return sparse.csr_matrix(mat / norms)


if __name__ == "__main__":
    from sys import argv
    main(argv)
