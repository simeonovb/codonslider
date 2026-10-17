"""Command line interface: codonslider max|gene --codon AGA --max-window 10 (--fasta FILE | --source | --genes)"""
from __future__ import annotations

import argparse
import logging
import sys

from .codonslider import geneslider, maxslider
from .version import __version__

ANALYSES = {'max': maxslider, 'gene': geneslider}


def _parser():
    parser = argparse.ArgumentParser(
        prog='codonslider',
        description='Find stretches of a codon in coding sequences with a growing sliding window. '
                    'Writes the results as CSV.')
    parser.add_argument('--version', action='version', version=f'%(prog)s {__version__}')
    parser.add_argument('analysis', choices=ANALYSES,
                        help="'max': the densest window per window size across all sequences (maxslider); "
                             "'gene': the densest window of every sequence per window size (geneslider)")
    parser.add_argument('-c', '--codon', required=True, help='the codon to count, e.g. AGA')
    parser.add_argument('--min-window', type=int, default=1, help='smallest window size in codons (default: 1)')
    parser.add_argument('--max-window', type=int, required=True, help='largest window size in codons')
    parser.add_argument('--step', type=int, default=1, help='codons the window moves per slide (default: 1)')

    sequences = parser.add_argument_group('sequences (choose one)').add_mutually_exclusive_group(required=True)
    sequences.add_argument('-f', '--fasta', help='your own CCDS FASTA file')
    sequences.add_argument('--source', choices=['ncbi', 'ensembl'],
                           help='download the NCBI CCDS release or the Ensembl CDS of all transcripts')
    sequences.add_argument('--genes', nargs='+', metavar='GENE',
                           help='retrieve only these genes (symbols or Ensembl IDs) from Ensembl')
    parser.add_argument('--species', default='human', help="'human', 'mouse' or an Ensembl species name "
                                                           "(default: human)")
    parser.add_argument('--assembly', help='human assembly: GRCh38 (default), GRCh37 or T2T')
    parser.add_argument('--cache-dir', help='where downloads are kept (default: ~/.codonslider)')
    parser.add_argument('--no-cache', action='store_true', help='delete the downloaded sequences afterwards')

    parser.add_argument('-o', '--output', help='CSV file to write (default: standard output)')
    parser.add_argument('-q', '--quiet', action='store_true', help='only report warnings and errors')
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = _parser()
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.WARNING if args.quiet else logging.INFO, format='%(message)s',
                        stream=sys.stderr)

    try:
        results, _, _ = ANALYSES[args.analysis](
            path_to_file=args.fasta, source=args.source, genes=args.genes, species=args.species,
            assembly=args.assembly, cache_dir=args.cache_dir, cache=not args.no_cache, steps=args.step,
            codon_of_interest=args.codon, initial_window_size=args.min_window, max_window_size=args.max_window)
    except (ValueError, OSError) as error:
        # Wrong arguments, missing files and network failures; the message says what to do
        print(f'codonslider: error: {error}', file=sys.stderr)
        return 1

    results.to_csv(args.output if args.output else sys.stdout, index=False)
    return 0
