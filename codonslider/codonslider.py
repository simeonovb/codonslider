from Bio import SeqIO
from datetime import datetime
import pandas as pd
from .sources import input_file, gene_info

def slidwin(elements, current_window_size, step, codon, max_window_size, identifier):
    # Window size and step are given in codons, the sequence is sliced in nucleotides
    window_length = current_window_size * 3
    if len(elements) <= window_length:
        return {'Window_Size': current_window_size, 'Codon_Frequency': 0, 'Sequence': elements}

    max_window_size = min(max_window_size, len(elements) // 3)  # Ensure the maximum window size is within the sequence length
    max_codon_freq = 0  # Track the maximum codon frequency across all records for the current window size
    sequence_with_max_codon_freq = ''

    for a in range(0, len(elements) - window_length + 1, step * 3):
        win = elements[a:a + window_length]

        # Count the codon of interest in the current window
        codons = [win[b:b + 3] for b in range(0, len(win), 3)]
        codon_freq = codons.count(codon) / current_window_size

        # Update the maximum codon frequency information if needed
        if codon_freq >= max_codon_freq:
            max_codon_freq = codon_freq
            sequence_with_max_codon_freq = win

    return {'Window_Size': current_window_size, 'Codon_Frequency': max_codon_freq, 'Sequence': sequence_with_max_codon_freq, 'Gene': identifier}


def maxslider(path_to_file=None, steps=1, codon_of_interest=None, initial_window_size=1, max_window_size=None,
              *, source=None, species='human', genes=None, assembly=None, cache_dir=None, cache=True):
    """
    Reads a FASTA file, iterates through sequence records, and applies the
    sliding window analysis (slidwin) for a range of window sizes.

    Args:
        path_to_file (str): Path to your own input CCDS FASTA file. Leave it out to
            retrieve the sequences with source/genes instead.
        steps (int): How many codons the window moves per slide.
        codon_of_interest (str): The codon to count (e.g., 'ATG').
        initial_window_size (int): The starting size for the sliding window, in codons.
        max_window_size (int): The maximum size for the sliding window, in codons.
        source (str): 'ncbi' to download the NCBI CCDS release (human or mouse) or 'ensembl'
            to download the CDS of all transcripts of the current Ensembl release.
        species (str): 'human', 'mouse' or, for Ensembl, any Ensembl species name (e.g. 'danio_rerio').
        genes (list): Analyse only these genes (symbols or Ensembl IDs), retrieved from the Ensembl REST API.
        assembly (str): Human genome assembly: 'GRCh38' (default), 'GRCh37' or 'T2T' (T2T-CHM13v2.0).
        cache_dir (str): Where downloaded files are stored. Defaults to ~/.codonslider.
        cache (bool): Set to False to delete the downloaded sequences as soon as the analysis is finished.

    Returns:
        tuple: A tuple containing (DataFrame of results, list of window sizes, list of max frequencies).
    """

    if codon_of_interest is None or max_window_size is None:
        raise TypeError("codon_of_interest and max_window_size are required.")
    with input_file(path_to_file, source=source, species=species, genes=genes, assembly=assembly,
                    cache_dir=cache_dir, cache=cache) as path_to_file:
        return _maxslider(path_to_file, steps, codon_of_interest, initial_window_size, max_window_size)


def _maxslider(path_to_file, steps, codon_of_interest, initial_window_size, max_window_size):
    # Initialize variables to keep track of the highest density information
    all_results = []
    current_window_size = initial_window_size
    window_sizes_for_plot = []
    max_codon_frequencies_for_plot = []

    # storing the current time in the variable
    date = datetime.now()
    # Displays Time
    current_time = date.strftime('%H:%M:%S')
    print("Starting time:", date)
    # Loop over different window sizes
    while current_window_size <= max_window_size:
        # NOTE: Opening the file inside the loop is necessary since SeqIO.parse
        # consumes the generator. You must open/close on each iteration.
        with open(path_to_file, mode='r') as handle:
            max_codon_freq_across_records = 0
            sequence_with_max_codon_freq_across_records = ''

            for record in SeqIO.parse(handle, 'fasta'):
                identifier = record.id
                sequence = record.seq

                # Skip records with sequences shorter than the current window size
                if len(sequence) < current_window_size * 3:
                    # print(f"Skipping record {identifier} with sequence length {len(sequence)}")
                    continue

                # Use the imported slidwin function
                active_info = slidwin(str(sequence), current_window_size, steps, codon_of_interest,
                                      max_window_size, identifier)

                # Update the maximum codon frequency information across all records
                if active_info['Codon_Frequency'] >= max_codon_freq_across_records:
                    max_codon_freq_across_records = active_info['Codon_Frequency']
                    sequence_with_max_codon_freq_across_records = active_info['Sequence']

            # Append the maximum codon frequency across all records for the current window size
            all_results.append({'Window_Size': current_window_size,
                                'Codon_Frequency': max_codon_freq_across_records,
                                'Sequence': sequence_with_max_codon_freq_across_records})

            # Append data for plotting
            window_sizes_for_plot.append(current_window_size)
            max_codon_frequencies_for_plot.append(max_codon_freq_across_records)

            print(
                f"Highest codon frequency across all records for window size {current_window_size} / Frequency: {max_codon_freq_across_records} / Gene: {identifier} / finished at: {datetime.now().strftime('%H:%M:%S')}")

        # Increase the window size by one codon
        current_window_size += 1

    final_result_df = pd.DataFrame(all_results)

    # Return the results needed by the remaining script
    return final_result_df, window_sizes_for_plot, max_codon_frequencies_for_plot


def geneslider(path_to_file=None, steps=1, codon_of_interest=None, initial_window_size=1, max_window_size=None,
               *, source=None, species='human', genes=None, assembly=None, cache_dir=None, cache=True):
    """
    Reads a FASTA file, iterates through sequence records, and applies the
    sliding_growing_window analysis for a range of window sizes.

    Args:
        path_to_file (str): Path to your own input CCDS FASTA file. Leave it out to
            retrieve the sequences with source/genes instead.
        steps (int): How many codons the window moves per slide.
        codon_of_interest (str): The codon to count (e.g., 'ATG').
        initial_window_size (int): The starting size for the sliding window, in codons.
        max_window_size (int): The maximum size for the sliding window, in codons.
        source (str): 'ncbi' to download the NCBI CCDS release (human or mouse) or 'ensembl'
            to download the CDS of all transcripts of the current Ensembl release.
        species (str): 'human', 'mouse' or, for Ensembl, any Ensembl species name (e.g. 'danio_rerio').
        genes (list): Analyse only these genes (symbols or Ensembl IDs), retrieved from the Ensembl REST API.
        assembly (str): Human genome assembly: 'GRCh38' (default), 'GRCh37' or 'T2T' (T2T-CHM13v2.0).
        cache_dir (str): Where downloaded files are stored. Defaults to ~/.codonslider.
        cache (bool): Set to False to delete the downloaded sequences as soon as the analysis is finished.

    Returns:
        tuple: A tuple containing (DataFrame of results, list of window sizes, list of max frequencies).
    """

    if codon_of_interest is None or max_window_size is None:
        raise TypeError("codon_of_interest and max_window_size are required.")
    with input_file(path_to_file, source=source, species=species, genes=genes, assembly=assembly,
                    cache_dir=cache_dir, cache=cache) as path_to_file:
        return _geneslider(path_to_file, steps, codon_of_interest, initial_window_size, max_window_size)


def _geneslider(path_to_file, steps, codon_of_interest, initial_window_size, max_window_size):
    # Initialize variables to keep track of the highest density information
    all_results = []
    current_window_size = initial_window_size
    window_sizes_for_plot = []
    max_codon_frequencies_for_plot = []

    # storing the current time in the variable
    date = datetime.now()
    # Displays Time
    current_time = date.strftime('%H:%M:%S')
    print("Starting time:", date)

    # Loop over different window sizes
    while current_window_size <= max_window_size:
        # Variables for tracking highest frequency in the current window size iteration
        frequencies_in_current_window = []

        # NOTE: The file must be opened inside the loop because SeqIO.parse consumes the generator
        with open(path_to_file, mode='r') as handle:
            for record in SeqIO.parse(handle, 'fasta'):
                identifier = record.id
                sequence = str(record.seq)  # Ensure sequence is a string for slicing

                # Skip sequences shorter than the current window size
                if len(sequence) < current_window_size * 3:
                    # print(f"Skipping protein {identifier} with sequence length {len(sequence)}")
                    continue

                # Apply the sliding window function
                active_info = slidwin(sequence, current_window_size, steps, codon_of_interest,
                                                     max_window_size, identifier)

                # Add the gene symbol and gene ID, if the FASTA header names the gene
                active_info['Gene_Symbol'], active_info['Gene_ID'] = gene_info(record.description)

                # Append result for the current sequence/window size combination
                all_results.append(active_info)

                # Keep track of the frequencies just for the current window size iteration
                frequencies_in_current_window.append(active_info['Codon_Frequency'])

        # --- End of iteration for the current_window_size ---

        # Update plotting lists with the highest frequency found in this window size
        if frequencies_in_current_window:
            max_freq_for_current_window = max(frequencies_in_current_window)
            window_sizes_for_plot.append(current_window_size)
            max_codon_frequencies_for_plot.append(max_freq_for_current_window)

            # Print information for the current window size
            print(
                f"Highest codon frequency across all proteins for window size {current_window_size} / Frequency: {max_freq_for_current_window} / Gene: {identifier} / finished at: {datetime.now().strftime('%H:%M:%S')}")
        else:
            print(f"No records processed for window size {current_window_size}.")

        # Increase the window size by one codon
        current_window_size += 1

    # Create the final DataFrame
    final_result_df = pd.DataFrame(all_results)

    # Return the results needed by the remaining script (for saving/plotting)
    return final_result_df, window_sizes_for_plot, max_codon_frequencies_for_plot
