"""
Example: where are the longest stretches of the AGA codon in the human CCDS?

Downloads the NCBI CCDS release (about 60 MB, kept in ~/.codonslider for the next run), saves the results
as CSV and plots the highest density per window size. Needs matplotlib: pip install matplotlib
"""
import logging

import matplotlib.pyplot as plt

from codonslider import geneslider, maxslider

# Show the progress of the downloads and the analysis
logging.basicConfig(level=logging.INFO, format='%(message)s')

codon = 'AGA'
sequences = {'source': 'ncbi', 'species': 'human'}  # or {'path_to_file': 'CCDS.fasta'} for your own file

# The densest window of each size across the whole transcriptome
densest, window_sizes, frequencies = maxslider(codon_of_interest=codon, initial_window_size=1, max_window_size=20,
                                               **sequences)
densest.to_csv('max_codon_stretches.csv', index=False)

# The best window of every gene, to see which genes carry the longest stretches
per_gene, _, _ = geneslider(codon_of_interest=codon, initial_window_size=10, max_window_size=10, **sequences)
per_gene.to_csv('codon_stretches_per_gene.csv', index=False)
print(per_gene.sort_values('Codon_Frequency', ascending=False).head(10))

plt.plot(window_sizes, frequencies, marker='o')
plt.title(codon)
plt.xlabel('Window Size (codons)')
plt.ylabel('Highest Codon Density')
plt.savefig('max_codon_stretches.png')
plt.show()
