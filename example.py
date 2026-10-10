import matplotlib.pyplot as plt
from codonslider import maxslider
from codonslider import geneslider

# Inputs
codon_of_interest = str(input("Which codon should be processed? ").upper())
initial_window_size = int(input("Initial size of the sliding window (in codons)? "))
max_window_size = int(input("Maximum size of the sliding window (in codons)? "))

# Analysis: 'maxslider' gives the highest codon density per window size across all sequences,
# 'geneslider' lists every gene with its best stretch per window size
analysis = input("Highest density per window size (1) or the stretches of every gene (2)? ").strip()
if analysis == '2':
    slider = geneslider
    output_file = "codon_stretches_per_gene.csv"
else:
    slider = maxslider
    output_file = "max_codon_stretches.csv"

# Sequences: either your own CCDS file (option 1) or a download from NCBI / Ensembl (option 2)
choice = input("Use your own CCDS file (1) or download the sequences (2)? ").strip()
if choice == '1':
    path_to_file = input("Path to your CCDS FASTA file? ").strip()
    source = None
    species = 'human'
    assembly = None
    cache = True
else:
    path_to_file = None
    source = input("Download from 'ncbi' (CCDS, human or mouse) or 'ensembl' (all CDS, any species)? ").strip().lower()
    species = input("Which species (e.g. human, mouse)? ").strip().lower()
    # The genome assembly can only be chosen for human
    assembly = None
    if species == 'human':
        assembly = input("Which assembly: GRCh38, GRCh37 or T2T? ").strip() or None
    # Keeping the file makes the next run faster; otherwise it is deleted when the analysis is finished
    cache = input("Keep the downloaded sequences for the next run (y/n)? ").strip().lower() != 'n'

# Run the codon slider. Both functions take the same arguments
final_result_df, window_sizes_for_plot, max_codon_frequencies_for_plot = slider(path_to_file=path_to_file,
                                                                                codon_of_interest=codon_of_interest,
                                                                                initial_window_size=initial_window_size,
                                                                                max_window_size=max_window_size,
                                                                                source=source, species=species, assembly=assembly,
                                                                                cache=cache)

# Save the final DataFrame to a CSV file
final_result_df.to_csv(output_file, index=False)

# With geneslider, show the 10 genes with the highest codon density in the largest window
if analysis == '2':
    largest_window = final_result_df[final_result_df["Window_Size"] == max_window_size]
    print(largest_window.sort_values("Codon_Frequency", ascending=False).head(10))

# Display and save the graph
plt.plot(window_sizes_for_plot, max_codon_frequencies_for_plot, marker='o')
plt.title(codon_of_interest)
plt.xlabel('Window Size (codons)')
plt.ylabel('Highest Codon Density')
plt.grid(False)
plt.savefig("Highest codon frequency across all records.tif", format="tif")
plt.show()