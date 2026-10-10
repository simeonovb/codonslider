# Codonslider -> codon stretch identificator

codonslider is a new package implemented in python to determine the percentage of codon stretches within certain region of the gene.
It works by applying the sliding window principle which scans the CCDS sequence and records the percentage of the codon of interest.
The sliding window can grow overtime and record the codon accumulation per window size. 
It is a useful tool for people who study protein translation and tRNA biology.

## Python Support

Python >=3.9 is required.

## Dependencies

1. Biopython
2. pandas

## Installation
**Using pip**

```python
pip install codonslider
```
## Contribution Guidelines

We welcome contributions to the software!

For reporting bugs or making suggestions, the most effective method is to raise an issue on the GitHub issue tracker. 
GitHub allows you to classify your issues, letting us know whether it's a bug report, feature request, or feedback for the authors.

If you would like to contribute changes to the code, please submit a [pull request](https://github.com/simeonovb/codonslider/pulls). 
For guidance on creating a pull request, please refer to the [documentation on pull requests](https://help.github.com/en/articles/about-pull-requests).

## Usage

codonslider has two functions. They take the same arguments and differ in what they report:

- `maxslider` reports, for each window size, the single window with the highest frequency of the codon across all sequences.
  Use it to see how long the densest stretch in the transcriptome is.
- `geneslider` reports, for each window size, the best window of every gene.
  Use it to find out which genes carry the stretches.

### maxslider: the highest codon density per window size

```python
from codonslider import maxslider
import matplotlib.pyplot as plt

path_to_file = "CCDS.fasta"
steps = 1 # the window moves one codon at a time
codon_of_interest = "AGA"
initial_window_size = 1 # in codons, i.e. 3 nucleotides
max_window_size = 4 # in codons, i.e. 12 nucleotides

final_result_df, window_sizes_for_plot, max_codon_frequencies_for_plot = maxslider(
                                                path_to_file = path_to_file,
                                                steps = steps, codon_of_interest = codon_of_interest,
                                                initial_window_size = initial_window_size,
                                                max_window_size = max_window_size)
# Save the final DataFrame to a CSV file
final_result_df.to_csv("codon_stretches.csv", index=False)

# Display and save the graph
plt.plot(window_sizes_for_plot, max_codon_frequencies_for_plot, marker='o')
plt.title(codon_of_interest)
plt.xlabel('Window Size (codons)')
plt.ylabel('Highest Codon Density')
plt.grid(False)
plt.savefig("codon_stretches.tif", format="tif")
plt.show()
```

The dataframe has one row per window size:

```
   Window_Size  Codon_Frequency   Sequence
0            1         1.000000        AGA
1            2         1.000000     AGAAGA
2            3         0.666667  AGAAGAGCC
```

### geneslider: the best stretch of every gene

```python
from codonslider import geneslider

final_result_df, window_sizes_for_plot, max_codon_frequencies_for_plot = geneslider(
                                                              path_to_file = "CCDS.fasta",
                                                              codon_of_interest = "AGA",
                                                              initial_window_size = 1,
                                                              max_window_size = 3)
# Save the stretches of all genes to a CSV file
final_result_df.to_csv("codon_stretches_per_gene.csv", index=False)

# The 10 genes with the highest density of the codon in a window of 3 codons
window = final_result_df[final_result_df["Window_Size"] == 3]
print(window.sort_values("Codon_Frequency", ascending=False).head(10))
```

The dataframe has one row per gene and window size, with three additional columns: `Gene` (the name of the
sequence in the FASTA file), `Gene_Symbol` and `Gene_ID`. For three genes (TP53, BRCA2 and MYC) it looks like this:

```
   Window_Size  Codon_Frequency   Sequence                   Gene Gene_Symbol          Gene_ID
0            1         1.000000        AGA   TP53|ENST00000269305        TP53  ENSG00000141510
1            1         1.000000        AGA  BRCA2|ENST00000380152       BRCA2  ENSG00000139618
2            1         1.000000        AGA    MYC|ENST00000621592         MYC  ENSG00000136997
3            2         0.500000     AGAGAC   TP53|ENST00000269305        TP53  ENSG00000141510
4            2         1.000000     AGAAGA  BRCA2|ENST00000380152       BRCA2  ENSG00000139618
5            2         0.500000     AGACAG    MYC|ENST00000621592         MYC  ENSG00000136997
6            3         0.333333  AGAGACCGG   TP53|ENST00000269305        TP53  ENSG00000141510
7            3         0.666667  AGAAGAGCC  BRCA2|ENST00000380152       BRCA2  ENSG00000139618
8            3         0.333333  AGACAGATC    MYC|ENST00000621592         MYC  ENSG00000136997
```

The gene symbol and gene ID are available for all sequences retrieved from NCBI or Ensembl. `Gene_ID` is the
identifier of the source: the NCBI Gene ID (e.g. `7157`) for NCBI and the Ensembl gene ID (e.g. `ENSG00000141510`)
for Ensembl. With your own file, the two columns are filled in when the FASTA headers name the gene in the Ensembl
(`gene:ENSG00000141510 gene_symbol:TP53`) or NCBI RefSeq (`[gene=TP53] [db_xref=GeneID:7157]`) format, and are
left empty otherwise.

Genes shorter than the window are left out for that window size. On a whole transcriptome the table gets large
(about 35,000 rows per window size for the human CCDS), so keep the range of window sizes small.

### Choosing the sequences

There are two ways to provide the sequences. They work the same for `maxslider` and `geneslider`.

**Option 1: your own CCDS file**

```python
maxslider(path_to_file = "CCDS.fasta", codon_of_interest = "AGA", initial_window_size = 1, max_window_size = 4)
```

**Option 2: retrieve the sequences from NCBI or Ensembl**

```python
# NCBI CCDS release (human or mouse)
maxslider(source = "ncbi", species = "human", codon_of_interest = "AGA", initial_window_size = 1, max_window_size = 4)

# CDS of all transcripts in the current Ensembl release (any Ensembl species, e.g. "danio_rerio")
maxslider(source = "ensembl", species = "mouse", codon_of_interest = "AGA", initial_window_size = 1, max_window_size = 4)

# Only selected genes (symbols or Ensembl IDs), retrieved through the Ensembl REST API
geneslider(genes = ["TP53", "BRCA2"], species = "human", codon_of_interest = "AGA", initial_window_size = 1, max_window_size = 4)
```

For human, the genome assembly can be chosen with `assembly`: `"GRCh38"` (default), `"GRCh37"` or `"T2T"` (T2T-CHM13v2.0).

```python
maxslider(source = "ncbi", assembly = "GRCh37", codon_of_interest = "AGA", initial_window_size = 1, max_window_size = 4)
maxslider(source = "ensembl", assembly = "T2T", codon_of_interest = "AGA", initial_window_size = 1, max_window_size = 4)
```

| | `source = "ncbi"` | `source = "ensembl"` | `genes = [...]` |
|---|---|---|---|
| `"GRCh38"` | CCDS, current release | CDS of all transcripts, current release | Ensembl REST API |
| `"GRCh37"` | CCDS, last GRCh37 release (Hs105, 2013) | CDS of all transcripts, Ensembl GRCh37 archive | Ensembl GRCh37 REST API |
| `"T2T"` | CDS of the RefSeq annotation (there is no CCDS for T2T) | CDS of all transcripts, Ensembl 2022_07 annotation | not available |

If you do not want to keep the downloaded files, use `cache = False`: the sequences are downloaded for this run only and deleted as soon as the analysis is finished.

```python
maxslider(source = "ncbi", cache = False, codon_of_interest = "AGA", initial_window_size = 1, max_window_size = 4)
```

To only download a file, or to refresh it, use `fetch_ccds` / `fetch_genes`, which return the path to the FASTA file:

```python
from codonslider import fetch_ccds, fetch_genes

path_to_file = fetch_ccds(source = "ncbi", species = "human", force = True) # force = True downloads it again
path_to_file = fetch_genes(["TP53", "BRCA2"], species = "human", output_file = "my_genes.fasta")
```

Note that the Ensembl file is not restricted to CCDS: it contains the CDS of every transcript, including
alternative isoforms and incomplete CDS, and is much larger (about 550 MB for human, versus 60 MB for NCBI CCDS).
For selected genes, the canonical transcript of each gene is used.

Where:

```
path_to_file: a CCDS fasta file which contains all sequences that needs to be analyzed (option 1)
source = "ncbi" or "ensembl", to download the sequences instead of using your own file (option 2)
species = "human", "mouse" or, for Ensembl, any Ensembl species name. Default is "human"
genes = a list of gene symbols or Ensembl IDs to analyse instead of the whole transcriptome
assembly = "GRCh38", "GRCh37" or "T2T". Human only. Default is "GRCh38"
cache = False deletes the downloaded sequences when the analysis is finished. Default is True (keep them)
steps = How many codons each slide moves. Default is 1
codon_of_interest = a codon you would like to analyse
initial_window_size = the initial size of the sliding window in codons (1 codon = 3 nucleotides). Minimum is 1
max_window_size = the final size of the growing sliding window in codons
```

Returns:

```
final_result_df = a dataframe with the window size (in codons), the highest codon frequency and the sequence of that window.
                  maxslider: one row per window size. geneslider: one row per gene and window size, with the gene symbol and gene ID
window_sizes_for_plot = the progression of the window expansion
max_codon_frequencies_for_plot = maximum codon accumulation recorded for the current window size (across all genes)
```

