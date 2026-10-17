"""Retrieval of coding sequences from NCBI and Ensembl."""
from __future__ import annotations

import gzip
import json
import logging
import os
import re
import shutil
import tempfile
import time
import urllib.error
import urllib.parse
import urllib.request
from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path

logger = logging.getLogger(__name__)

NCBI_CCDS_URL = 'https://ftp.ncbi.nlm.nih.gov/pub/CCDS/{release}/CCDS_nucleotide.current.fna.gz'
# Last CCDS release built on GRCh37 (GRCh37.p13, October 2013)
NCBI_GRCH37_RELEASE = 'archive/Hs105'
# There is no CCDS for T2T-CHM13, so the CDS of the RefSeq annotation is used
NCBI_T2T_URL = ('https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/'
                'GCF_009914755.1_T2T-CHM13v2.0_cds_from_genomic.fna.gz')
ENSEMBL_CDS_URL = 'https://ftp.ensembl.org/pub/current_fasta/{species}/cds/'
ENSEMBL_GRCH37_CDS_URL = 'https://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/cds/'
ENSEMBL_T2T_URL = ('https://ftp.ensembl.org/pub/rapid-release/species/Homo_sapiens/GCA_009914755.4/ensembl/geneset/'
                   '2022_07/Homo_sapiens-GCA_009914755.4-2022_07-cds.fa.gz')
# Tables used to add the gene to FASTA headers that do not name it
NCBI_CCDS_TABLE_URL = 'https://ftp.ncbi.nlm.nih.gov/pub/CCDS/{release}/CCDS.current.txt'
ENSEMBL_GRCH37_SYMBOLS_URL = ('https://grch37.ensembl.org/biomart/martservice?query=' + urllib.parse.quote(
    '<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query virtualSchemaName="default" formatter="TSV" '
    'header="0" uniqueRows="1" completionStamp="1"><Dataset name="hsapiens_gene_ensembl" interface="default">'
    '<Attribute name="ensembl_gene_id"/><Attribute name="external_gene_name"/></Dataset></Query>'))
ENSEMBL_REST_URL = 'https://rest.ensembl.org'
ENSEMBL_GRCH37_REST_URL = 'https://grch37.rest.ensembl.org'

# NCBI only publishes CCDS for these two species
NCBI_SPECIES = {'human': 'human', 'homo_sapiens': 'human', 'mouse': 'mouse', 'mus_musculus': 'mouse'}
# Common names accepted for Ensembl; any other Ensembl species name (e.g. 'danio_rerio') is passed through
ENSEMBL_SPECIES = {'human': 'homo_sapiens', 'mouse': 'mus_musculus'}
# Accepted names for the human genome assemblies
ASSEMBLIES = {'grch38': 'GRCh38', 'hg38': 'GRCh38', 'grch37': 'GRCh37', 'hg19': 'GRCh37',
              't2t': 'T2T', 't2t-chm13': 'T2T', 't2t-chm13v2.0': 'T2T', 'chm13': 'T2T', 'hs1': 'T2T'}

DEFAULT_CACHE_DIR = Path.home() / '.codonslider'
# Ensembl REST limit for POST sequence/id
REST_BATCH_SIZE = 50


def _open_url(url, data=None):
    headers = {'User-Agent': 'codonslider'}
    if data is not None:
        data = json.dumps(data).encode()
        headers.update({'Content-Type': 'application/json', 'Accept': 'application/json'})
    request = urllib.request.Request(url, data=data, headers=headers)
    try:
        return urllib.request.urlopen(request, timeout=60)
    except urllib.error.HTTPError as error:
        if error.code == 404:
            raise FileNotFoundError(f"{url} does not exist") from error
        raise ConnectionError(f"Request to {url} failed: {error.code} {error.reason}") from error
    except urllib.error.URLError as error:
        raise ConnectionError(f"Could not reach {url}: {error.reason}") from error


def _download_gz(url, destination, annotate=None):
    """
    Downloads a gzipped FASTA file and stores it decompressed at destination.
    annotate is an optional function which returns the text to add to a header line.
    """
    logger.info("Downloading %s", url)
    # Write to a temporary file first, so an interrupted download never leaves a truncated FASTA in the cache
    partial = destination.with_name(destination.name + '.part')
    try:
        with _open_url(url) as response, gzip.open(response) as unzipped, open(partial, mode='wb') as handle:
            if annotate is None:
                shutil.copyfileobj(unzipped, handle)
            else:
                for line in unzipped:
                    if line.startswith(b'>'):
                        line = line.rstrip() + annotate(line.decode()).encode() + b'\n'
                    handle.write(line)
        os.replace(partial, destination)
    finally:
        partial.unlink(missing_ok=True)
    logger.info("Saved to %s", destination)


def _ccds_annotation(release):
    """The CCDS FASTA headers only hold the CCDS ID, so the gene is looked up in the CCDS table of the release."""
    genes = {}
    with _open_url(NCBI_CCDS_TABLE_URL.format(release=release)) as response:
        for line in response:
            if not line.startswith(b'#'):
                columns = line.decode().split('\t')
                genes[columns[4]] = f' gene:{columns[3]} gene_symbol:{columns[2]}'
    # Header: >CCDS2.2|Hs110|chr1
    return lambda header: genes.get(header[1:].split('|')[0], '')


def _ensembl_grch37_annotation():
    """The Ensembl GRCh37 FASTA headers hold the gene ID but not the symbol, which is looked up in BioMart."""
    symbols = {}
    # BioMart regularly answers with an error page or an incomplete table, so check its completion stamp and retry
    for _attempt in range(3):
        try:
            with _open_url(ENSEMBL_GRCH37_SYMBOLS_URL) as response:
                lines = response.read().decode().splitlines()
        except OSError:
            lines = []
        if lines and lines[-1] == '[success]':
            for line in lines[:-1]:
                gene_id, _, symbol = line.partition('\t')
                if symbol:
                    symbols[gene_id] = f' gene_symbol:{symbol}'
            break
        time.sleep(5)
    else:
        logger.warning("Could not retrieve the gene symbols for Ensembl GRCh37 (Ensembl BioMart did not answer), so "
                       "only the gene IDs are included. Download again later with fetch_ccds(..., force=True) to "
                       "add the symbols.")

    def annotate(header):
        gene_id = re.search(r'\bgene:([^.\s]+)', header)
        return symbols.get(gene_id.group(1), '') if gene_id else ''
    return annotate


def gene_info(description: str) -> tuple[str, str]:
    """
    Reads the gene symbol and gene ID from the header of a FASTA record. Understands the headers of
    Ensembl ('gene:ENSG00000141510.18 gene_symbol:TP53') and of NCBI RefSeq ('[gene=TP53] [db_xref=GeneID:7157]').

    Returns:
        tuple: (gene symbol, gene ID). An empty string is used for what the header does not contain.
    """
    symbol = re.search(r'\bgene_symbol:(\S+)|\[gene=([^\]]+)\]', description)
    gene_id = re.search(r'\bgene:(\S+)|\bGeneID:(\d+)', description)
    return (symbol.group(1) or symbol.group(2) if symbol else '',
            gene_id.group(1) or gene_id.group(2) if gene_id else '')


def _ensembl_species(species):
    species = species.strip().lower().replace(' ', '_')
    return ENSEMBL_SPECIES.get(species, species)


def _human_assembly(assembly, species):
    """Returns 'GRCh37', 'T2T' or None, which stands for the current assembly (GRCh38 for human)."""
    if assembly is None:
        return None
    name = ASSEMBLIES.get(assembly.strip().lower())
    if name is None:
        raise ValueError(f"Unknown assembly '{assembly}'. Use 'GRCh38', 'GRCh37' or 'T2T'.")
    if _ensembl_species(species) != 'homo_sapiens':
        raise ValueError(f"The assembly can only be chosen for human, not '{species}'.")
    return None if name == 'GRCh38' else name


def _ensembl_post(endpoint, payload, assembly=None):
    rest_url = ENSEMBL_GRCH37_REST_URL if assembly == 'GRCh37' else ENSEMBL_REST_URL
    with _open_url(rest_url + endpoint, data=payload) as response:
        return json.load(response)


def fetch_ccds(source: str = 'ncbi', species: str = 'human', assembly: str | None = None,
               cache_dir: str | os.PathLike | None = None, force: bool = False) -> str:
    """
    Downloads the coding sequences of a species and returns the path to the FASTA file.
    The file is kept in the cache directory and reused on the next call.

    Args:
        source (str): 'ncbi' for the NCBI CCDS release (human or mouse only) or
            'ensembl' for the CDS of all transcripts in the current Ensembl release.
        species (str): 'human', 'mouse' or, for Ensembl, any Ensembl species name (e.g. 'danio_rerio').
        assembly (str): Human genome assembly: 'GRCh38' (default), 'GRCh37' or 'T2T' (T2T-CHM13v2.0).
            NCBI has no CCDS for T2T, so the CDS of its RefSeq annotation is downloaded instead.
        cache_dir (str): Where the downloaded file is stored. Defaults to ~/.codonslider.
        force (bool): Download again even if the file is already in the cache.

    Returns:
        str: Path to the downloaded FASTA file.
    """
    source = source.lower()
    assembly = _human_assembly(assembly, species)
    cache_dir = Path(cache_dir) if cache_dir else DEFAULT_CACHE_DIR
    cache_dir.mkdir(parents=True, exist_ok=True)

    if source == 'ncbi':
        ncbi_species = NCBI_SPECIES.get(species.strip().lower().replace(' ', '_'))
        if ncbi_species is None:
            raise ValueError(f"NCBI CCDS is only available for human and mouse, not '{species}'. "
                             f"Use source='ensembl' for other species.")
        if assembly == 'T2T':
            logger.info("NCBI has no CCDS release for T2T-CHM13. Using the CDS of the NCBI RefSeq annotation instead.")
            url = NCBI_T2T_URL
            destination = cache_dir / 'ncbi_human_T2T_RefSeq_cds.fasta'
        elif assembly == 'GRCh37':
            release = NCBI_GRCH37_RELEASE
            destination = cache_dir / 'ncbi_human_GRCh37_CCDS_nucleotide.fasta'
        else:
            release = f'current_{ncbi_species}'
            destination = cache_dir / f'ncbi_{ncbi_species}_CCDS_nucleotide.current.fasta'
        if force or not destination.exists():
            if assembly == 'T2T':
                _download_gz(url, destination)
            else:
                _download_gz(NCBI_CCDS_URL.format(release=release), destination, annotate=_ccds_annotation(release))

    elif source == 'ensembl':
        ensembl_species = _ensembl_species(species)
        if assembly == 'T2T':
            destination = cache_dir / 'ensembl_homo_sapiens_T2T_cds.all.fasta'
            if force or not destination.exists():
                _download_gz(ENSEMBL_T2T_URL, destination)
        else:
            if assembly == 'GRCh37':
                directory = ENSEMBL_GRCH37_CDS_URL
                destination = cache_dir / 'ensembl_homo_sapiens_GRCh37_cds.all.fasta'
            else:
                directory = ENSEMBL_CDS_URL.format(species=ensembl_species)
                destination = cache_dir / f'ensembl_{ensembl_species}_cds.all.fasta'
            if force or not destination.exists():
                # The file name contains the assembly (e.g. Homo_sapiens.GRCh38.cds.all.fa.gz), so look it up
                try:
                    with _open_url(directory) as response:
                        listing = response.read().decode()
                except FileNotFoundError as error:
                    raise ValueError(f"No Ensembl CDS file found for species '{species}'. "
                                     f"Use the Ensembl species name, e.g. 'homo_sapiens'.") from error
                file_name = re.search(r'href="([^"]+\.cds\.all\.fa\.gz)"', listing)
                if file_name is None:
                    raise ValueError(f"No Ensembl CDS file found for species '{species}'.")
                _download_gz(directory + file_name.group(1), destination,
                             annotate=_ensembl_grch37_annotation() if assembly == 'GRCh37' else None)

    else:
        raise ValueError(f"Unknown source '{source}'. Use 'ncbi' or 'ensembl'.")

    return str(destination)


def fetch_genes(genes: str | list[str], species: str = 'human', assembly: str | None = None,
                output_file: str | os.PathLike | None = None, cache_dir: str | os.PathLike | None = None) -> str:
    """
    Retrieves the CDS of the canonical transcript of the given genes through the
    Ensembl REST API and returns the path to a FASTA file containing them.

    Args:
        genes (list): Gene symbols (e.g. 'TP53') or Ensembl gene/transcript IDs.
        species (str): 'human', 'mouse' or any Ensembl species name. Only used for gene symbols.
        assembly (str): Human genome assembly: 'GRCh38' (default) or 'GRCh37'. T2T is not available here.
        output_file (str): Where to write the FASTA file. Defaults to a file in the cache directory.
        cache_dir (str): Directory used when output_file is not given. Defaults to ~/.codonslider.

    Returns:
        str: Path to the FASTA file.
    """
    if isinstance(genes, str):
        genes = [genes]
    genes = list(dict.fromkeys(genes))
    if not genes:
        raise ValueError("No genes given.")
    ensembl_species = _ensembl_species(species)
    assembly = _human_assembly(assembly, species)
    if assembly == 'T2T':
        raise ValueError("The Ensembl REST API does not offer T2T. Individual genes can be retrieved for "
                         "GRCh38 or GRCh37; for T2T, download all sequences with source='ncbi' or 'ensembl'.")

    # Resolve each gene to a transcript ID: first as a symbol, then as an Ensembl ID
    transcripts = {}
    found = _ensembl_post(f'/lookup/symbol/{ensembl_species}', {'symbols': genes}, assembly)
    remaining = [gene for gene in genes if gene not in found]
    if remaining:
        by_id = _ensembl_post('/lookup/id', {'ids': remaining}, assembly)
        found.update({gene: entry for gene, entry in by_id.items() if entry})
    for gene, entry in found.items():
        name = entry.get('display_name') or gene
        if entry.get('object_type') == 'Transcript':
            # A transcript is named after its gene (BRCA1-203)
            transcript, gene_id, symbol = entry['id'], entry.get('Parent', ''), name.rsplit('-', 1)[0]
        else:
            transcript, gene_id, symbol = entry.get('canonical_transcript'), entry['id'], name
        if transcript:
            # Remove the version number (ENST00000269305.9 -> ENST00000269305)
            transcripts[transcript.split('.')[0]] = (name, f'gene:{gene_id} gene_symbol:{symbol}')

    not_found = [gene for gene in genes if gene not in found]
    if not_found:
        raise ValueError(f"Not found in Ensembl ({ensembl_species}): {', '.join(not_found)}")

    transcript_ids = list(transcripts)
    records = []
    for a in range(0, len(transcript_ids), REST_BATCH_SIZE):
        try:
            batch = transcript_ids[a:a + REST_BATCH_SIZE]
            records.extend(_ensembl_post('/sequence/id?type=cds', {'ids': batch}, assembly))
        except FileNotFoundError:
            # Ensembl answers 404 when none of the transcripts in the batch has a CDS (non-coding genes)
            continue
    if not records:
        raise ValueError(f"None of the given genes has a coding sequence: {', '.join(genes)}")
    missing = set(transcript_ids) - {record['query'] for record in records}
    if missing:
        logger.warning("Skipping genes without a coding sequence: %s",
                       ', '.join(transcripts[m][0] for m in sorted(missing)))

    if output_file is None:
        cache_dir = Path(cache_dir) if cache_dir else DEFAULT_CACHE_DIR
        cache_dir.mkdir(parents=True, exist_ok=True)
        output_file = cache_dir / f"ensembl_{ensembl_species}{'_' + assembly if assembly else ''}_selected_genes.fasta"
    with open(output_file, mode='w') as handle:
        for record in records:
            name, gene = transcripts[record['query']]
            handle.write(f">{name}|{record['id']} {gene}\n{record['seq']}\n")

    return str(output_file)


def resolve_input(path_to_file: str | None = None, source: str | None = None, species: str = 'human',
                  genes: list[str] | None = None, assembly: str | None = None,
                  cache_dir: str | os.PathLike | None = None) -> str:
    """Returns the FASTA file to analyse: the user's own file, or one retrieved from NCBI/Ensembl."""
    if path_to_file is not None:
        if source is not None or genes is not None or assembly is not None:
            raise ValueError("Give either path_to_file or source/genes/assembly, not both.")
        return path_to_file
    if genes is not None:
        if source not in (None, 'ensembl'):
            raise ValueError("Individual genes can only be retrieved from Ensembl (source='ensembl').")
        return fetch_genes(genes, species=species, assembly=assembly, cache_dir=cache_dir)
    if source is None:
        raise ValueError("Give path_to_file (your own CCDS FASTA file) or source ('ncbi' or 'ensembl').")
    return fetch_ccds(source=source, species=species, assembly=assembly, cache_dir=cache_dir)


@contextmanager
def input_file(path_to_file: str | None = None, source: str | None = None, species: str = 'human',
               genes: list[str] | None = None, assembly: str | None = None,
               cache_dir: str | os.PathLike | None = None, cache: bool = True) -> Iterator[str]:
    """Provides the FASTA file to analyse. With cache=False, retrieved sequences are deleted afterwards."""
    if cache or path_to_file is not None:
        yield resolve_input(path_to_file, source=source, species=species, genes=genes, assembly=assembly,
                            cache_dir=cache_dir)
        return
    if cache_dir is not None:
        raise ValueError("cache_dir cannot be combined with cache=False.")
    # Download into a temporary directory, which is removed even if the analysis fails
    with tempfile.TemporaryDirectory(prefix='codonslider_') as temporary_dir:
        yield resolve_input(source=source, species=species, genes=genes, assembly=assembly, cache_dir=temporary_dir)
    logger.info("Deleted the downloaded sequences.")
