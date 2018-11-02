


# Semantic Discovery via Singular Value Decomposition

This is a summer project done with the Illinois Geometry Lab at the University
of Illinois.

## Project Description

This project attempted to apply [Singular Value Decomposition][SVD] to gene
sequence data to perform some analysis.

[SVD]: https://en.wikipedia.org/wiki/Singular_value_decomposition


## How to Build and Run

### Prerequisites

	* [Python 3][Py3]: The program is written in Python.
	* [Numpy and Scipy][NPSP]: High-performance matrix data type and linear
		algebra routines making use of those matrices
	* [Biopython][BP]: A library for working with genome data
	* [Matplotlib][MP]: Used to create plots and graphs

[Py3]: https://python.org
[NPSP]: https://scipy.org
[BP]: https://biopython.org
[MP]: https://matplotlib.org

### Running the Program

At a command prompt, run

```bash
$ python3 main.py
```

to download genome data (large), construct the term-document matrix, and create plots.


## Notes about the program

### Data Sources

To operate, the program requires genome files in the GenBank format. There are
two ways to obtain those files.

	* Local directories
	* According to a `genomes.txt` file

Local directories are recursively searched to find `.gb` and `.gbk` files. A
directory containing a file called `genomes.txt` is considered a *remote directory*.
The `genomes.txt` file (see below for details) references files on servers at
the National Center for Biotechnology Information, that will be retrieved over
FTP if not already present.

Right now, to specify the directory you will need to edit the call to
`get_matrix` made in `main`.

Currently, the main driver of the program (the stuff related to finding and
acquiring the data files) is rather messy. The biggest symptom of this is the
existence of both `get_matrix()` and `get_matrix2()`. You should probably prefer
using `get_matrix2()` as it is more cleanly written, has simpler control flow,
and less edge cases.

### The `genomes.txt` Format

The format of a `genomes.txt` file is quite simple.

	* Blank lines and lines where the first non-whitespace character is `#` are
		ignored.
	* The first significant line contains a number of `key=value` pairs, separated
		by whitespace. Currently, the only meaningful attribute is `destdir`, which
		specifies the output directory of the `genomes.txt` file.
	* Genomes are specified over two lines:

		1. The path to the genome on the NCBI FTP server
		2. The name of the file to where the genome will be saved.

If the server-path ends with `.gz`, the file will be decompressed after downloading.

