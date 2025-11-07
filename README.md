# Strainberry 2

Strainberry2 is the new and improved version of Strainberry. It is conceived to perform strain-level metagenome assembly with modern long reads (e.g., PacBio HiFi, ONT R10.4).


## Installation

Strainberry has been developed and tested under a Linux environment.

### System prerequisites

- Rust compiler >= 1.82

### Build from source

Clone the repository:

```
git clone https://github.com/rvicedomini/strainberry2.git
```

To build and install the `strainberry2` executable in your system:

```
cd strainberry2
cargo install --path .
```

### Software dependencies

Strainberry2 requires the following bioinformatics tool to be available in the system:
- samtools >= 1.15
- minimap2 >= 2.27

A simple way to install these dependencies would be to create a conda environment `strainberry2` as follows:
```
conda create -n strainberry2 -c bioconda 'samtools>=1.15' 'minimap2>=2.27'
```

Then to run strainberry2:
```
conda run -n strainberry2 strainberry2 [OPTIONS] --reference <PATH> --out-dir <PATH>
```

## Quick usage

Strainberry2 requires two input files:
- **A set of long reads**, preferably in FASTQ format as base qualities are used.
- **A reference assembly**, in either FASTA or GFA format. Consider using a strain-oblivious assembler such as metaMDBG or metaFlye to generate such input. 

All input files could be provided either uncompressed or gzipped.

To run Strainberry2 with **ONT R10.4** reads using 8 threads:
```
strainberry2 --in-ont reads.fastq.gz -r assembly.fasta --out-dir output -t 8
```

To run Strainberry2 with **PacBio HiFi** reads using 8 threads:
```
strainberry2 --in-hifi reads.fastq.gz -r assembly.fasta --out-dir output -t 8
```

### Other useful parameters

The following is a list of parameters that you might want to tune.
The full list is available running `strainberry2 --help`, which also shows all default values.

- `-q <NUM>`: minimum MAPQ value to consider a read alignment
- `--min-snv <NUM>`: minimum number of phased SNVs to retain a haplotype
- `--min-alt-count <NUM>`: minimum number of alternative-allele observations
- `--min-alt-frac <FLOAT>`: minimum fraction of alternative-allele observations

## Output files

The strain-aware output assembly is written in the `assembly.fasta` file in inside the output directory provided with the `--out-dir` parameter. A corresponding assembly graph in GFA format is also available in the `assembly.gfa` file.

## Example

In order to verify that Strainberry2 has been correctly installed, it is possible to test it on a tiny dataset in the `example` sub-directory.

Given a strain-oblivious assembly (e.g., file `assembly.fa.gz`) and a collection of accurate ONT reads (e.g., file `reads.fa.gz`), it is possible to run Strainberry2 using 4 threads as follows:

```
cd example
strainberry2 --in-ont reads.fa.gz -r assembly.fa.gz -o output -t 4
```

If everything worked as expected, you will find the strain-aware assembly `assembly.fasta` and the corresponding strain-aware assembly graph `assembly.gfa` in the `output` directory. Both files should contain 5 contigs. Each one of them corresponds to the 5 *Escherichia coli* strain sequences present in the input data and collapsed in the input single-contig assembly. 

## Citation

A preprint will be available soon.

