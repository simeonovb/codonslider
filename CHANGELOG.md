# Changelog

## 3.0.0 (unreleased)

### Fixed

- A sequence exactly as long as the window was given a codon frequency of 0. It is now analysed like any other
  sequence (with geneslider, its row also lacked the gene name).
- The progress messages named the last gene in the file instead of the gene with the densest window.
- maxslider crashed on an empty FASTA file. It now returns an empty table.
- The codon was not matched in lower case sequences, and a lower case codon (`'aga'`) matched nothing.
- A sequence without the codon in any window was reported with the last window as its `Sequence` (often the stop
  codon, e.g. `TAA` for ISG15 when searching for AGA). `Sequence` is now empty when the frequency is 0, and in
  maxslider so are the gene columns, since no gene has the codon in a window of that size.

### Changed

- maxslider and geneslider are about 60 to 200 times faster: the FASTA file is read once instead of once per window
  size, and the codons in a window are counted from a cumulative count instead of one by one. The results are the
  same as before, apart from the fixes above.
- maxslider also reports the gene with the densest window: new columns `Gene`, `Gene_Symbol` and `Gene_ID`.
- maxslider leaves out window sizes longer than every sequence, as geneslider does, instead of reporting them with a
  frequency of 0 and an empty sequence.
- Progress and download messages go through the `logging` module instead of `print`, so they can be turned off or
  redirected. Use `logging.basicConfig(level=logging.INFO)` to see them.
- **Breaking:** `slidwin(sequence, window_size, step, codon, identifier=None)` no longer takes the unused
  `max_window_size` argument. A sequence shorter than the window now gives an empty `Sequence`.
- The arguments are checked: the codon must be three of the letters A, C, G and T (or U), and the window sizes and
  step must be whole numbers of at least 1, with `max_window_size` not smaller than `initial_window_size`.
- The package is built from `pyproject.toml` (`setup.py` is removed), declares NumPy as a dependency, and includes
  type hints.

### Added

- The `codonslider` command line tool (also `python -m codonslider`).
- Tests (100% coverage, without network access) and continuous integration on GitHub Actions, with linting by ruff.

## 2.0 and earlier

See the [commit history](https://github.com/simeonovb/codonslider/commits/master).
