"""Sliding window search for stretches of one codon in coding sequences."""
from __future__ import annotations

import logging
from collections import defaultdict
from collections.abc import Iterable

import numpy as np
import pandas as pd
from Bio import SeqIO

from .sources import gene_info, input_file

logger = logging.getLogger(__name__)

# maxslider reports no gene: many genes often share the highest frequency, and geneslider lists them
MAXSLIDER_COLUMNS = ['Window_Size', 'Codon_Frequency', 'Sequence']
GENESLIDER_COLUMNS = MAXSLIDER_COLUMNS + ['Gene', 'Gene_Symbol', 'Gene_ID']
NUCLEOTIDES = set('ACGT')


def _check_codon(codon):
    """Returns the codon in upper case DNA letters ('aga' and 'AGA' both give 'AGA', 'AGU' gives 'AGT')."""
    if not isinstance(codon, str):
        raise TypeError(f"The codon must be a string such as 'AGA', not {codon!r}.")
    normalized = codon.strip().upper().replace('U', 'T')
    if len(normalized) != 3 or not set(normalized) <= NUCLEOTIDES:
        raise ValueError(f"'{codon}' is not a codon. Use three of the letters A, C, G and T (or U).")
    return normalized


def _check_positive(name, value):
    if isinstance(value, bool) or not isinstance(value, (int, np.integer)) or value < 1:
        raise ValueError(f"{name} must be a whole number of at least 1, not {value!r}.")


def _codon_counts(sequence, codon):
    """
    Cumulative count of the codon in the reading frame of the sequence: counts[i] is the number of times
    it occurs in the first i codons, so a window of codons [a, b) contains counts[b] - counts[a] of them.
    """
    n_codons = len(sequence) // 3
    codons = np.frombuffer(sequence[:n_codons * 3].encode('ascii'), dtype='S3')
    counts = np.zeros(n_codons + 1, dtype=np.int64)
    np.cumsum(codons == codon.encode('ascii'), out=counts[1:])
    return counts


def _densest_window(sequence, counts, window_size, step):
    """
    Returns (frequency, window sequence) of the densest window, the last one if several are equally dense.
    A sequence shorter than the window, or without the codon in any window, gives (0.0, '').
    """
    n_codons = len(counts) - 1
    if n_codons < window_size:
        return 0.0, ''
    # Number of codons of interest in the windows starting at codon 0, step, 2 * step, ...
    per_window = counts[window_size::step] - counts[:n_codons - window_size + 1:step]
    best = len(per_window) - 1 - int(np.argmax(per_window[::-1]))
    if per_window[best] == 0:
        # Every window ties at 0, so no window is worth showing
        return 0.0, ''
    start = best * step * 3
    return int(per_window[best]) / window_size, sequence[start:start + window_size * 3]


def slidwin(sequence: str, window_size: int, step: int, codon: str, identifier: str | None = None) -> dict:
    """
    Slides a window over one coding sequence and finds the window with the highest frequency of the codon.
    Only codons in the reading frame of the sequence are counted.

    Args:
        sequence (str): The coding sequence, starting with its first codon.
        window_size (int): The size of the window, in codons.
        step (int): How many codons the window moves per slide.
        codon (str): The codon to count (e.g. 'AGA').
        identifier (str): Name of the sequence, copied to the result.

    Returns:
        dict: Window_Size, Codon_Frequency (the fraction of the codons in the window that are the codon of
        interest), Sequence (the densest window, the last one if several are equally dense) and Gene (identifier).
        A sequence shorter than the window, or without the codon in any window, has a frequency of 0 and an
        empty Sequence.
    """
    codon = _check_codon(codon)
    _check_positive('window_size', window_size)
    _check_positive('step', step)
    sequence = sequence.upper()
    frequency, window = _densest_window(sequence, _codon_counts(sequence, codon), window_size, step)
    return {'Window_Size': window_size, 'Codon_Frequency': frequency, 'Sequence': window, 'Gene': identifier}


def _scan(path_to_file, codon, step, window_sizes):
    """
    Reads the FASTA file once and yields the densest window of every sequence for every window size.
    Window sizes larger than a sequence are skipped for that sequence.
    """
    with open(path_to_file) as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            sequence = str(record.seq).upper()
            counts = _codon_counts(sequence, codon)
            symbol, gene_id = gene_info(record.description)
            for window_size in window_sizes:
                if window_size > len(counts) - 1:
                    break
                frequency, window = _densest_window(sequence, counts, window_size, step)
                yield {'Window_Size': window_size, 'Codon_Frequency': frequency, 'Sequence': window,
                       'Gene': record.id, 'Gene_Symbol': symbol, 'Gene_ID': gene_id}


def _densest_per_window_size(windows: Iterable[dict]) -> dict[int, list[dict]]:
    """maxslider: keeps the densest window across all sequences per window size (the last one if tied)."""
    densest: dict[int, dict] = {}
    for window in windows:
        current = densest.get(window['Window_Size'])
        if current is None or window['Codon_Frequency'] >= current['Codon_Frequency']:
            densest[window['Window_Size']] = window
    return {window_size: [window] for window_size, window in densest.items()}


def _all_per_window_size(windows: Iterable[dict]) -> dict[int, list[dict]]:
    """geneslider: keeps the densest window of every sequence per window size."""
    by_window_size = defaultdict(list)
    for window in windows:
        by_window_size[window['Window_Size']].append(window)
    return by_window_size


def _run(collect, columns, path_to_file, steps, codon_of_interest, initial_window_size, max_window_size,
         **retrieval):
    if codon_of_interest is None or max_window_size is None:
        raise TypeError("codon_of_interest and max_window_size are required.")
    codon = _check_codon(codon_of_interest)
    _check_positive('steps', steps)
    _check_positive('initial_window_size', initial_window_size)
    _check_positive('max_window_size', max_window_size)
    if max_window_size < initial_window_size:
        raise ValueError(f"max_window_size ({max_window_size}) is smaller than "
                         f"initial_window_size ({initial_window_size}).")
    window_sizes = range(initial_window_size, max_window_size + 1)

    with input_file(path_to_file, **retrieval) as path_to_file:
        logger.info("Searching %s for stretches of %s in windows of %d to %d codons",
                    path_to_file, codon, initial_window_size, max_window_size)
        by_window_size = collect(_scan(path_to_file, codon, steps, window_sizes))

    rows, window_sizes_for_plot, max_codon_frequencies_for_plot = [], [], []
    for window_size in window_sizes:
        windows = by_window_size.get(window_size)
        if not windows:
            logger.info("Window size %d: no sequence is that long", window_size)
            continue
        highest = max(window['Codon_Frequency'] for window in windows)
        if highest:
            logger.info("Window size %d: highest frequency %.3f", window_size, highest)
        else:
            logger.info("Window size %d: no window contains %s", window_size, codon)
        rows.extend(windows)
        window_sizes_for_plot.append(window_size)
        max_codon_frequencies_for_plot.append(highest)

    # Columns not listed (the gene of maxslider's window) are left out
    return pd.DataFrame(rows, columns=columns), window_sizes_for_plot, max_codon_frequencies_for_plot


def maxslider(path_to_file: str | None = None, steps: int = 1, codon_of_interest: str | None = None,
              initial_window_size: int = 1, max_window_size: int | None = None, *, source: str | None = None,
              species: str = 'human', genes: list[str] | None = None, assembly: str | None = None,
              cache_dir: str | None = None, cache: bool = True) -> tuple[pd.DataFrame, list[int], list[float]]:
    """
    Finds, for each window size, the highest frequency of the codon in any window across all sequences:
    the densest possible stretch. Use geneslider to see which genes reach it.

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
        tuple: (DataFrame with one row per window size, list of window sizes, list of max frequencies).
        The DataFrame has the columns Window_Size, Codon_Frequency and Sequence (one of the densest windows,
        the last one in the file; empty if no window contains the codon). Window sizes larger than every
        sequence are left out.
    """
    return _run(_densest_per_window_size, MAXSLIDER_COLUMNS, path_to_file, steps, codon_of_interest,
                initial_window_size, max_window_size, source=source, species=species, genes=genes,
                assembly=assembly, cache_dir=cache_dir, cache=cache)


def geneslider(path_to_file: str | None = None, steps: int = 1, codon_of_interest: str | None = None,
               initial_window_size: int = 1, max_window_size: int | None = None, *, source: str | None = None,
               species: str = 'human', genes: list[str] | None = None, assembly: str | None = None,
               cache_dir: str | None = None, cache: bool = True) -> tuple[pd.DataFrame, list[int], list[float]]:
    """
    Finds, for each window size, the window with the highest frequency of the codon in every sequence.

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
        tuple: (DataFrame with one row per sequence and window size, list of window sizes,
        list of max frequencies). Sequences shorter than the window are left out for that window size.
    """
    return _run(_all_per_window_size, GENESLIDER_COLUMNS, path_to_file, steps, codon_of_interest,
                initial_window_size, max_window_size, source=source, species=species, genes=genes,
                assembly=assembly, cache_dir=cache_dir, cache=cache)
