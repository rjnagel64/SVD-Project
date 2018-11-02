
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
#   (done?)

# TODO: Create visualization of the rank 2 approximation using matplotlib. (done?)

# TODO: One of the sequences in the drosophila melanogaster genome (sequence 0)
# has a length that is not a multiple of three. Biopython suggests appending an
# 'N' to the end of the sequence to rectify this.

# TODO: Replace ad-hoc `genomes.txt` format with JSON/YAML? Eh. Not worth it.
#       (Also, for some reason, the python standard library does not come with
#       a module for YAML.)
# TODO: Allow genomes.txt to specify its output directory, so genomes.txt can be
#       checked into version control, but the genomes are not. Alternatively,
#       add a .gitignore to the directory containing the genomes.txt file.
#       - The first non-blank line should contain key=value pairs, split on
#         whitespace.

# TODO: Proper argument parsing
#   * Extract the the logic from `get_matrix()`
#   * `-l`, `--local` [directory]
#   * `-r`, `--remote` [genomes.txt]

# TODO: In doc comments, I talk about `(mat, proteins, genomes)` at least 3 times.
#       It should become its own data type, probably. (So it has a name.)

DNA_SEQUENCE_DIR = "dna-sequences"
PROTEIN_SEQUENCE_DIR = "protein-sequences"
MATRIX_SAVE_DIR = "matrix"
PLOTS_DIR = "plots"
DEFAULT_GENOME_PATHS_FILE = "genomes.txt"

def main(args):
    # (mat, proteins, genomes) = get_matrix(args)
    (mat, proteins, genomes) = get_matrix2("./new-dna/", is_remote=True)

    print(mat.shape)

    (W, ss, Vh) = sparse.linalg.svds(mat, k=2, which="LM")
    print(ss)
    print(Vh)

    mat = mat
    xr = (-25, 25)
    yr = (-10, 300)
    create_plot(mat, ss, Vh, name="not_normalized.pdf", xrange=xr, yrange=yr)

    mat2 = inverse_frequency(mat)
    xr = (-50, 50)
    yr = (-10, 500)
    create_plot(mat2, ss, Vh, name="tf_idf_normalized.pdf", xrange=xr, yrange=yr)

    mat3 = norm(mat)
    xr = (-2, 2)
    yr = (-2, 2)
    create_plot(mat3, ss, Vh, name="row_normalized.eps", xrange=xr, yrange=yr)

# TODO: Remove argument `ss` from `create_plot`. It is unused.
def create_plot(mat, ss, Vh, name, xrange=None, yrange=None):
    """Create a scatter plot of words in the new basis, and save it to a PDF file

    arguments:
    - mat: The term-document matrix
    - ss: Singular values (not used? may remove this parameter)
    - Vh: The matrix on the right in SVD. It forms a basis for V
    - name: The name of the plot. The result will be saved to `PLOTS_DIR/{name}`.
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

    plot_path = os.path.join(PLOTS_DIR, name)
    print(f"Saving plot to {plot_path}...")
    plt.savefig(plot_path)

# TODO: Move command-line argument parsing elsewhere, and separate local/remote data collection.
def get_matrix(args):
    """Obtain the term-document matrix, either by loading NCBI FTP paths from a
    file, all (.gbk) paths in a subdirectory, or creating the matrix from
    scratch.
    """
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
            get_ncbi_data(DEFAULT_GENOME_PATHS_FILE)

        if not os.path.exists(PROTEIN_SEQUENCE_DIR):
            print(f"Writing proteins...")
            paths = list(map(lambda name: os.path.join(DNA_SEQUENCE_DIR, name),
                             os.listdir(DNA_SEQUENCE_DIR)))
            write_proteins(paths)

    print("Constructing matrix...")
    (mat, proteins, genomes) = create_matrix2()

    print(f"Matrix created: {mat.shape}.")

    print("Saving data...")
    save_data(mat, proteins, genomes)

    return (mat, proteins, genomes)

def get_matrix2(path, is_remote):
    """Get the term-document matrix.

    If `MATRIX_SAVE_DIR` exists, the matrix will be loaded from there.
    Otherwise, the matrix will be created and saved to `MATRIX_SAVE_DIR`.

    arguments:
    - path: A file or directory path. Its interpretation depends on the value of
        `is_remote`.
    - is_remote: If `is_remote` is `False`, `path` is a path to a directory
        containing `.gb` and `.gbk` files. It will be recursively searched.
        If `is_remote` is `True`, `path` is the path to a `genomes.txt` file
        specifying genome files on NCBI's GenBank FTP server, which will be
        retrieved if they are not already present.

    returns:
    TODO: Actually write out what the return type of this function is. (see `load_data`)
    - The term-document matrix, genome labels, and protein labels
    """
    # If the matrix already exists, just load it.
    if os.path.exists(MATRIX_SAVE_DIR):
        return load_data(MATRIX_SAVE_DIR)

    # Otherwise, the matrix must be constructed.
    print("Constructing matrix...")

    # Search for the sequence files that will be used to construct the matrix.
    print("Searching for GenBank files...")
    if is_remote:
        dir = load_remote(path)
    else:
        dir = path

    # Find the genome files to use.
    paths = find_files(dir)
    print(f"Found {len(paths)} GenBank files.")

    # For each genome, determine if the protein sequences have already been
    # extracted. If not, write the proteins.
    to_write = []
    for path in paths:
        (_, filename) = os.path.split(path)
        (base, _ext) = os.path.splitext(filename)

        protein_path = os.path.join(PROTEIN_SEQUENCE_DIR, base)
        if not os.path.exists(protein_path):
            to_write.append(path)

    if to_write != []:
        print(f"Writing {len(to_write)} protein sequence files...")
        write_proteins(to_write)

    else:
        print("No protein sequences to write.")

    # Now that all the protein sequences have been created, construct the matrix.
    (mat, proteins, genomes) = create_matrix2()

    # Save the matrix.
    save_data(mat, proteins, genomes)

    return (mat, proteins, genomes)

def save_data(mat, proteins, genomes):
    """Save the term-document matrix into the subdirectory `MATRIX_SAVE_DIR`.

    arguments:
    - mat: A scipy sparse array.
    - proteins: A list of protein sequences such that the `i`th row
        of `mat` corresponds to `proteins[i]`.
    - genomes: A list of genome names such that the `j`th column of `mat`
        corresponds to `genomes[j]`.
    """
    if not os.path.exists(MATRIX_SAVE_DIR):
        os.mkdir(MATRIX_SAVE_DIR)

    #  np.save(os.path.join(MATRIX_SAVE_DIR, "matrix.npy"), mat)
    sparse.save_npz(os.path.join(MATRIX_SAVE_DIR, "matrix.npz"), mat)

    with open(os.path.join(MATRIX_SAVE_DIR, "proteins.txt"), "w") as f:
        f.writelines(proteins)

    with open(os.path.join(MATRIX_SAVE_DIR, "genomes.txt"), "w") as f:
        f.writelines(genomes)

def load_data(dir):
    """Load the term-document matrix from a directory.

    arguments:
    - dir: A directory path containing files `matrix.npz`, `proteins.txt`, and
      `genomes.txt`

    returns:
    A tuple `(mat, proteins, genomes)` such that:
    - `mat` is a scipy sparse array.
    - `proteins` is the list of all protein sequences such that the `i`th row
        of `mat` corresponds to `proteins[i]`.
    - `genomes` is the list of genome names such that the `j`th column of `mat`
        corresponds to `genomes[j]`.
    """
    print(f"Loading matrix from save...")
    mat = sparse.load_npz(os.path.join(dir, "matrix.npz"))

    print(f"Loading proteins from save...")
    with open(os.path.join(dir, "proteins.txt"), "r") as f:
        proteins = list(f)

    print(f"Loading genomes from save...")
    #Note that the filename used here is not to be confused with the default
    #file listing FTP paths on the NCBI server
    with open(os.path.join(dir, "genomes.txt"), "r") as f:
        genomes = list(f)

    print(f"Matrix loaded.")
    return (mat, proteins, genomes)

def load_remote(genomes_file):
    """Use a `genomes.txt` file to populate a directory with genomes.

    Files mentioned in `genomes.txt` that are not already present in the
    directory will be downloaded via FTP from NCBI's GenBank. The files
    will be saved in directory specified the 'dest_dir' attribute of the
    genomes file (if it exists), or the same directory as the genomes file
    otherwise. This directory will be created if it does not already exist.

    arguments:
    - genomes_file: The path to a `genomes.txt` file.

    returns:
    - The path of directory containing the downloaded files.
    """
    #  genomes_file = os.path.join(dir, "genomes.txt")
    try:
        with open(genomes_file, "r") as f:
            genomes = f.read()
    except FileNotFoundError:
        print(f"Error: Could not open {genomes_file}")
        raise RuntimeError()

    # Take only the nonempty, noncomment lines of the file.
    lines = list(filter(lambda l: l != "" and l[0] != "#",
                        map(str.strip,
                            genomes.splitlines())))

    # The first line contains key=value pairs.
    attrs_line = lines[0]
    attrs = dict(map(lambda pair: pair.split("="), attrs_line.split()))
    # This turns ['a', 'b', 'c', 'd', ...] into [('a', 'b'), ('c', 'd'), ...].
    files = list(zip(lines[1::2], lines[2::2]))

    # Find the destination directory, creating it if necessary.
    if attrs["dest_dir"]:
        dest_dir = attrs["dest_dir"]
    else:
        dest_dir = os.path.dirname(genomes_file)

    if not os.path.exists(dest_dir):
        os.mkdir(dest_dir)

    # Only download files that are not already present in `dest_dir`.
    to_retrieve = []
    for (path, file) in files:
        if not os.path.exists(os.path.join(dest_dir, file)):
            to_retrieve.append((path, file))

    print(f"{genomes_file} contains {len(files)} files, of which " +
          f"{len(to_retrieve)} must be downloaded.")

    # Retrieve the files and store them in the destination directory.
    retrieve_files(to_retrieve, dest_dir)

    return dest_dir

def find_files(dir):
    """Recursively search a directory for `.gb` and `.gbk` files.
    
    arguments:
    - dir: The directory to search.
    
    returns:
    A list of paths to files.
    """
    paths = []
    if os.path.isdir(dir):
        print(f"Checking {dir} for .gb and .gbk files...")

        for (path, _dirnames, filenames) in os.walk(dir):
            for filename in filenames:
                (name, ext) = os.path.splitext(filename)
                if ext in (".gb", ".gbk"):
                    paths.append(os.path.join(path, filename))

    else:
        # Aargh. Error reporting.
        print(f"Error: {dir} was supposed to be a directory.")
        raise RuntimeError()

    return paths

def get_ncbi_data(path):
    """Retrieve data from the NCBI servers and databanks.

    The data files retrieved are stored in `DNA_SEQUENCE_DIR`.

    arguments:
    - path: the path to the `genomes.txt` file
    """

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

    # Some files retrieved from the server may be gzip-compressed.
    to_decompress = []

    # Connect to the GenBank server.
    f = FTP("ftp.ncbi.nlm.nih.gov")
    f.login()

    # Retrieve the files.
    for (fp, fn) in filepaths:
        (path, filename) = os.path.split(fp)

        (_, ext) = os.path.splitext(filename)
        if ext == ".gz":
            to_decompress.append((filename, fn))
            dest_file = filename
        else:
            dest_file = fn

        print(f"Downloading {filename} ...")
        # Return to the root directory of the available files.
        f.cwd("/")
        # Go to the directory containing the file we want.
        f.cwd(path)

        with open(os.path.join(dest_dir, dest_file), "wb") as fd:
            f.retrbinary("RETR " + filename, fd.write)

    f.quit()

    print("Files downloaded.")

    # Decompress downloaded gzip files.
    for (cfn, dfn) in to_decompress:
        print(f"Decompressing {cfn}...")
        compressed_filename = os.path.join(dest_dir, cfn)
        decompressed_filename = os.path.join(dest_dir, dfn)

        # Decompress the file and copy its contents into a new file.
        with gzip.open(compressed_filename, "rb") as compressed_file:
            with open(decompressed_filename, "wb") as decompressed_file:
                copyfileobj(compressed_file, decompressed_file)

        # Delete the compressed file. It's not useful.
        os.remove(compressed_filename)

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

def write_proteins(paths):
    """Write the proteins from a given list of paths to .gbk files into files
    in the protein sequence directory.

    If `PROTEIN_SEQUENCE_DIR` does not exist, it will be created.

    arguments:
    - paths: A list of paths to `.gb` or `.gbk` files.
    """
    if not os.path.exists(PROTEIN_SEQUENCE_DIR):
        os.mkdir(PROTEIN_SEQUENCE_DIR)

    totals = Counter()
    for path in paths:
        (_root, filename) = os.path.split(path)
        (basename, _ext) = os.path.splitext(filename)

        print(f"Processing file {path}...")

        counter = Counter()
        with open(os.path.join(PROTEIN_SEQUENCE_DIR, basename), "w") as f:
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


# TODO: We aren't actually using this function; it should probably be eliminated.
def create_matrix():
    # This function should probably return the matrix variable.
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
        corresponds to `genomes[j]`.
    """
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
    inverse document frequency of that word and document.
    """
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
    `mat` divided by that row's norm.
    """

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
