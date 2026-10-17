import gzip
import importlib
import io
import json
import os
import urllib.error

import pytest

from codonslider import fetch_ccds, fetch_genes, geneslider, maxslider, sources
from codonslider.sources import gene_info, input_file, resolve_input

# The module codonslider.codonslider (the package exports functions of the same name)
core = importlib.import_module('codonslider.codonslider')


@pytest.fixture(autouse=True)
def no_home_cache(monkeypatch, tmp_path):
    """Keeps every test away from the real ~/.codonslider."""
    monkeypatch.setattr(sources, 'DEFAULT_CACHE_DIR', tmp_path / 'default_cache')


@pytest.fixture
def fake_urls(monkeypatch):
    """
    Replaces the network with a dict of URL -> response bytes, and records the requested URLs.
    A list of responses is answered in turn, one per request.
    """
    responses = {}
    requested = []

    def open_url(url, data=None):
        requested.append(url)
        if url not in responses:
            raise FileNotFoundError(f"{url} does not exist")
        response = responses[url]
        if isinstance(response, list):
            response = response.pop(0)
        return io.BytesIO(response)

    monkeypatch.setattr(sources, '_open_url', open_url)
    monkeypatch.setattr(sources.time, 'sleep', lambda seconds: None)
    return responses, requested


@pytest.fixture
def fake_ensembl(monkeypatch):
    """
    Replaces the Ensembl REST API with a dict of endpoint -> function(payload) returning the JSON answer.
    The calls are recorded as (endpoint, payload, assembly) in endpoints['calls'].
    """
    calls = []
    endpoints = {'calls': calls}

    def ensembl_post(endpoint, payload, assembly=None):
        calls.append((endpoint, payload, assembly))
        if endpoint not in endpoints:
            raise FileNotFoundError(endpoint)
        return endpoints[endpoint](payload)

    monkeypatch.setattr(sources, '_ensembl_post', ensembl_post)
    return endpoints


def lookup_genes(payload):
    """Ensembl symbol lookup that knows every gene and gives it the transcript ENST_<symbol>."""
    return {symbol: {'id': f'ENSG_{symbol}', 'display_name': symbol, 'object_type': 'Gene',
                     'canonical_transcript': f'ENST_{symbol}.1'} for symbol in payload['symbols']}


def cds_of(payload, skip=()):
    """Ensembl sequence endpoint that returns a CDS for every transcript except those in skip."""
    return [{'query': query, 'id': query + '.1', 'seq': 'ATGAGA'} for query in payload['ids'] if query not in skip]


class TestOpenUrl:
    @pytest.fixture
    def urlopen(self, monkeypatch):
        """Replaces urllib's urlopen; set urlopen.error to make it raise. The requests are recorded."""
        def fake(request, timeout):
            fake.requests.append(request)
            if fake.error:
                raise fake.error
            return io.BytesIO(b'ok')
        fake.requests = []
        fake.error = None
        monkeypatch.setattr(sources.urllib.request, 'urlopen', fake)
        return fake

    def test_get(self, urlopen):
        assert sources._open_url('https://example.org/file').read() == b'ok'
        request = urlopen.requests[0]
        assert request.get_method() == 'GET'
        assert request.get_header('User-agent') == 'codonslider'

    def test_post_sends_json(self, urlopen):
        sources._open_url('https://example.org/api', data={'ids': ['A']})
        request = urlopen.requests[0]
        assert request.get_method() == 'POST'
        assert json.loads(request.data) == {'ids': ['A']}
        assert request.get_header('Content-type') == 'application/json'

    def test_404_is_file_not_found(self, urlopen):
        urlopen.error = urllib.error.HTTPError('https://example.org/file', 404, 'Not Found', {}, None)
        with pytest.raises(FileNotFoundError, match='does not exist'):
            sources._open_url('https://example.org/file')

    def test_other_http_error_is_connection_error(self, urlopen):
        urlopen.error = urllib.error.HTTPError('https://example.org/file', 503, 'Service Unavailable', {}, None)
        with pytest.raises(ConnectionError, match='503'):
            sources._open_url('https://example.org/file')

    def test_unreachable_is_connection_error(self, urlopen):
        urlopen.error = urllib.error.URLError('Name or service not known')
        with pytest.raises(ConnectionError, match='Could not reach'):
            sources._open_url('https://example.org/file')


class TestEnsemblPost:
    @pytest.fixture
    def posted(self, monkeypatch):
        posted = []

        def open_url(url, data=None):
            posted.append((url, data))
            return io.BytesIO(b'{"TP53": {"id": "ENSG00000141510"}}')
        monkeypatch.setattr(sources, '_open_url', open_url)
        return posted

    def test_current_assembly(self, posted):
        answer = sources._ensembl_post('/lookup/symbol/homo_sapiens', {'symbols': ['TP53']})
        assert answer == {'TP53': {'id': 'ENSG00000141510'}}
        assert posted == [('https://rest.ensembl.org/lookup/symbol/homo_sapiens', {'symbols': ['TP53']})]

    def test_grch37(self, posted):
        sources._ensembl_post('/lookup/id', {'ids': ['ENSG00000141510']}, 'GRCh37')
        assert posted[0][0] == 'https://grch37.rest.ensembl.org/lookup/id'


class TestGeneInfo:
    def test_ensembl_header(self):
        header = 'ENST00000269305.9 cds chromosome:GRCh38:17:1:2:-1 gene:ENSG00000141510.18 gene_symbol:TP53'
        assert gene_info(header) == ('TP53', 'ENSG00000141510.18')

    def test_ncbi_header(self):
        header = 'lcl|NC_060941.1_cds_1 [gene=TP53] [db_xref=GeneID:7157] [protein=cellular tumor antigen p53]'
        assert gene_info(header) == ('TP53', '7157')

    def test_unknown_header(self):
        assert gene_info('my_sequence') == ('', '')


class TestAssembly:
    @pytest.mark.parametrize('name, expected', [(None, None), ('GRCh38', None), ('hg38', None),
                                                ('grch37', 'GRCh37'), ('hg19', 'GRCh37'),
                                                ('T2T', 'T2T'), ('chm13', 'T2T')])
    def test_names(self, name, expected):
        assert sources._human_assembly(name, 'human') == expected

    def test_unknown_assembly(self):
        with pytest.raises(ValueError, match='Unknown assembly'):
            sources._human_assembly('GRCh36', 'human')

    def test_only_for_human(self):
        with pytest.raises(ValueError, match='only be chosen for human'):
            sources._human_assembly('GRCh37', 'mouse')

    @pytest.mark.parametrize('species, expected', [('human', 'homo_sapiens'), ('Mouse', 'mus_musculus'),
                                                   ('Danio rerio', 'danio_rerio')])
    def test_ensembl_species(self, species, expected):
        assert sources._ensembl_species(species) == expected


class TestFetchCcds:
    def test_ncbi_adds_gene_to_headers(self, fake_urls, tmp_path):
        responses, _ = fake_urls
        release = 'current_human'
        responses[sources.NCBI_CCDS_URL.format(release=release)] = gzip.compress(
            b'>CCDS2.2|Hs110|chr1\nATGAGA\n>CCDS3.1|Hs110|chr1\nATGCCC\n')
        responses[sources.NCBI_CCDS_TABLE_URL.format(release=release)] = (
            b'#chromosome\tnc_accession\tgene\tgene_id\tccds_id\tccds_status\tcds_strand\n'
            b'1\tNC_000001.11\tSAMD11\t148398\tCCDS2.2\tPublic\t+\n')

        path = fetch_ccds(source='ncbi', species='human', cache_dir=tmp_path)

        assert path == str(tmp_path / 'ncbi_human_CCDS_nucleotide.current.fasta')
        with open(path) as handle:
            assert handle.read() == ('>CCDS2.2|Hs110|chr1 gene:148398 gene_symbol:SAMD11\nATGAGA\n'
                                     '>CCDS3.1|Hs110|chr1\nATGCCC\n')
        assert not os.path.exists(path + '.part')

    def test_cached_file_is_reused(self, fake_urls, tmp_path):
        _, requested = fake_urls
        cached = tmp_path / 'ncbi_mouse_CCDS_nucleotide.current.fasta'
        cached.write_text('>CCDS1.1\nATG\n')

        assert fetch_ccds(source='ncbi', species='mouse', cache_dir=tmp_path) == str(cached)
        assert requested == []

    def test_force_downloads_again(self, fake_urls, tmp_path):
        _, requested = fake_urls
        (tmp_path / 'ncbi_mouse_CCDS_nucleotide.current.fasta').write_text('>CCDS1.1\nATG\n')
        with pytest.raises(FileNotFoundError):
            fetch_ccds(source='ncbi', species='mouse', cache_dir=tmp_path, force=True)
        assert requested

    def test_ensembl_finds_file_in_listing(self, fake_urls, tmp_path):
        responses, _ = fake_urls
        directory = sources.ENSEMBL_CDS_URL.format(species='danio_rerio')
        responses[directory] = b'<a href="Danio_rerio.GRCz11.cds.all.fa.gz">Danio_rerio.GRCz11.cds.all.fa.gz</a>'
        responses[directory + 'Danio_rerio.GRCz11.cds.all.fa.gz'] = gzip.compress(b'>ENSDART1 gene:ENSDARG1\nATG\n')

        path = fetch_ccds(source='ensembl', species='danio_rerio', cache_dir=tmp_path)

        assert path == str(tmp_path / 'ensembl_danio_rerio_cds.all.fasta')
        with open(path) as handle:
            assert handle.read() == '>ENSDART1 gene:ENSDARG1\nATG\n'

    def test_ensembl_unknown_species(self, fake_urls, tmp_path):
        with pytest.raises(ValueError, match='No Ensembl CDS file found'):
            fetch_ccds(source='ensembl', species='unicorn', cache_dir=tmp_path)

    def test_ensembl_listing_without_cds_file(self, fake_urls, tmp_path):
        responses, _ = fake_urls
        responses[sources.ENSEMBL_CDS_URL.format(species='danio_rerio')] = b'<a href="README">README</a>'
        with pytest.raises(ValueError, match='No Ensembl CDS file found'):
            fetch_ccds(source='ensembl', species='danio_rerio', cache_dir=tmp_path)

    def test_default_cache_dir(self, fake_urls):
        responses, _ = fake_urls
        directory = sources.ENSEMBL_CDS_URL.format(species='danio_rerio')
        responses[directory] = b'<a href="Danio_rerio.GRCz11.cds.all.fa.gz">'
        responses[directory + 'Danio_rerio.GRCz11.cds.all.fa.gz'] = gzip.compress(b'>ENSDART1\nATG\n')

        path = fetch_ccds(source='ensembl', species='danio_rerio')
        assert path == str(sources.DEFAULT_CACHE_DIR / 'ensembl_danio_rerio_cds.all.fasta')

    def test_failed_download_leaves_no_file(self, fake_urls, tmp_path):
        responses, _ = fake_urls
        directory = sources.ENSEMBL_CDS_URL.format(species='danio_rerio')
        responses[directory] = b'<a href="Danio_rerio.GRCz11.cds.all.fa.gz">'
        responses[directory + 'Danio_rerio.GRCz11.cds.all.fa.gz'] = b'not gzip'

        with pytest.raises(OSError):
            fetch_ccds(source='ensembl', species='danio_rerio', cache_dir=tmp_path)
        assert list(tmp_path.iterdir()) == []

    def test_ncbi_only_human_and_mouse(self, tmp_path):
        with pytest.raises(ValueError, match='only available for human and mouse'):
            fetch_ccds(source='ncbi', species='danio_rerio', cache_dir=tmp_path)

    def test_unknown_source(self, tmp_path):
        with pytest.raises(ValueError, match='Unknown source'):
            fetch_ccds(source='genbank', cache_dir=tmp_path)


class TestFetchCcdsAssemblies:
    def test_ncbi_grch37_uses_last_grch37_release(self, fake_urls, tmp_path):
        responses, requested = fake_urls
        release = sources.NCBI_GRCH37_RELEASE
        responses[sources.NCBI_CCDS_URL.format(release=release)] = gzip.compress(b'>CCDS2.2|Hs105|chr1\nATG\n')
        responses[sources.NCBI_CCDS_TABLE_URL.format(release=release)] = (
            b'1\tNC_000001.10\tSAMD11\t148398\tCCDS2.2\tPublic\t+\n')

        path = fetch_ccds(source='ncbi', assembly='hg19', cache_dir=tmp_path)

        assert path == str(tmp_path / 'ncbi_human_GRCh37_CCDS_nucleotide.fasta')
        with open(path) as handle:
            assert handle.readline() == '>CCDS2.2|Hs105|chr1 gene:148398 gene_symbol:SAMD11\n'
        assert all('archive/Hs105' in url for url in requested)

    def test_ncbi_t2t_uses_refseq_cds(self, fake_urls, tmp_path):
        responses, requested = fake_urls
        header = b'>lcl|NC_060941.1_cds_1 [gene=TP53] [db_xref=GeneID:7157]\nATG\n'
        responses[sources.NCBI_T2T_URL] = gzip.compress(header)

        path = fetch_ccds(source='ncbi', assembly='T2T', cache_dir=tmp_path)

        assert path == str(tmp_path / 'ncbi_human_T2T_RefSeq_cds.fasta')
        with open(path, mode='rb') as handle:
            assert handle.read() == header
        # The RefSeq headers already name the gene, so no CCDS table is downloaded
        assert requested == [sources.NCBI_T2T_URL]

    def test_ensembl_t2t(self, fake_urls, tmp_path):
        responses, _ = fake_urls
        responses[sources.ENSEMBL_T2T_URL] = gzip.compress(b'>ENST1 gene:ENSG1 gene_symbol:TP53\nATG\n')

        path = fetch_ccds(source='ensembl', assembly='T2T', cache_dir=tmp_path)

        assert path == str(tmp_path / 'ensembl_homo_sapiens_T2T_cds.all.fasta')

    def test_ensembl_t2t_cached(self, fake_urls, tmp_path):
        _, requested = fake_urls
        (tmp_path / 'ensembl_homo_sapiens_T2T_cds.all.fasta').write_text('>ENST1\nATG\n')
        fetch_ccds(source='ensembl', assembly='T2T', cache_dir=tmp_path)
        assert requested == []

    @pytest.fixture
    def grch37_files(self, fake_urls):
        responses, requested = fake_urls
        directory = sources.ENSEMBL_GRCH37_CDS_URL
        responses[directory] = b'<a href="Homo_sapiens.GRCh37.75.cds.all.fa.gz">'
        responses[directory + 'Homo_sapiens.GRCh37.75.cds.all.fa.gz'] = gzip.compress(
            b'>ENST00000269305 cds:known chromosome:GRCh37:17:1:2:-1 gene:ENSG00000141510.11\nATG\n'
            b'>ENST00000000001 cds:known chromosome:GRCh37:1:1:2:1 gene:ENSG00000000001\nATG\n')
        return responses, requested

    def test_ensembl_grch37_adds_symbols_from_biomart(self, grch37_files, tmp_path):
        responses, _ = grch37_files
        responses[sources.ENSEMBL_GRCH37_SYMBOLS_URL] = b'ENSG00000141510\tTP53\nENSG00000000001\t\n[success]\n'

        path = fetch_ccds(source='ensembl', assembly='GRCh37', cache_dir=tmp_path)

        assert path == str(tmp_path / 'ensembl_homo_sapiens_GRCh37_cds.all.fasta')
        with open(path) as handle:
            headers = [line for line in handle if line.startswith('>')]
        assert headers[0].endswith('gene:ENSG00000141510.11 gene_symbol:TP53\n')
        # A gene without a symbol in BioMart keeps its header unchanged
        assert headers[1].endswith('gene:ENSG00000000001\n')

    def test_biomart_is_retried(self, grch37_files, tmp_path):
        responses, requested = grch37_files
        # BioMart first answers with a table cut short (no completion stamp), then with the full table
        responses[sources.ENSEMBL_GRCH37_SYMBOLS_URL] = [b'ENSG00000141510\tTP53\n',
                                                         b'ENSG00000141510\tTP53\n[success]\n']

        path = fetch_ccds(source='ensembl', assembly='GRCh37', cache_dir=tmp_path)

        assert requested.count(sources.ENSEMBL_GRCH37_SYMBOLS_URL) == 2
        with open(path) as handle:
            assert 'gene_symbol:TP53' in handle.readline()

    def test_biomart_down_still_downloads(self, grch37_files, tmp_path, caplog):
        _, requested = grch37_files

        path = fetch_ccds(source='ensembl', assembly='GRCh37', cache_dir=tmp_path)

        assert requested.count(sources.ENSEMBL_GRCH37_SYMBOLS_URL) == 3
        assert 'Could not retrieve the gene symbols' in caplog.text
        with open(path) as handle:
            assert handle.readline().endswith('gene:ENSG00000141510.11\n')


class TestFetchGenes:
    def test_writes_canonical_cds(self, fake_ensembl, tmp_path):
        fake_ensembl['/lookup/symbol/homo_sapiens'] = lambda payload: {
            'TP53': {'id': 'ENSG00000141510', 'display_name': 'TP53', 'object_type': 'Gene',
                     'canonical_transcript': 'ENST00000269305.9'}}
        fake_ensembl['/lookup/id'] = lambda payload: {
            'ENST00000357654': {'id': 'ENST00000357654', 'display_name': 'BRCA1-203', 'object_type': 'Transcript',
                                'Parent': 'ENSG00000012048'}}
        fake_ensembl['/sequence/id?type=cds'] = lambda payload: [
            {'query': query, 'id': query + '.1', 'seq': 'ATGAGA'} for query in payload['ids']]

        output_file = tmp_path / 'genes.fasta'
        path = fetch_genes(['TP53', 'ENST00000357654'], output_file=output_file)

        assert path == str(output_file)
        assert output_file.read_text() == (
            '>TP53|ENST00000269305.1 gene:ENSG00000141510 gene_symbol:TP53\nATGAGA\n'
            '>BRCA1-203|ENST00000357654.1 gene:ENSG00000012048 gene_symbol:BRCA1\nATGAGA\n')

    def test_unknown_gene(self, fake_ensembl, tmp_path):
        fake_ensembl['/lookup/symbol/homo_sapiens'] = lambda payload: {}
        fake_ensembl['/lookup/id'] = lambda payload: {'NOTAGENE': None}
        with pytest.raises(ValueError, match='Not found in Ensembl.*NOTAGENE'):
            fetch_genes(['NOTAGENE'], cache_dir=tmp_path)

    def test_no_coding_sequence(self, fake_ensembl, tmp_path):
        fake_ensembl['/lookup/symbol/homo_sapiens'] = lambda payload: {
            'MALAT1': {'id': 'ENSG00000251562', 'object_type': 'Gene', 'canonical_transcript': 'ENST00000534336.1'}}
        with pytest.raises(ValueError, match='coding sequence'):
            fetch_genes('MALAT1', cache_dir=tmp_path)

    def test_requests_sequences_in_batches(self, fake_ensembl, tmp_path):
        genes = [f'GENE{n}' for n in range(sources.REST_BATCH_SIZE * 2 + 20)]
        fake_ensembl['/lookup/symbol/homo_sapiens'] = lookup_genes
        fake_ensembl['/sequence/id?type=cds'] = cds_of

        path = fetch_genes(genes, cache_dir=tmp_path)

        batches = [payload['ids'] for endpoint, payload, _ in fake_ensembl['calls'] if endpoint.startswith('/sequence')]
        assert [len(batch) for batch in batches] == [50, 50, 20]
        with open(path) as handle:
            assert sum(line.startswith('>') for line in handle) == len(genes)

    def test_skips_genes_without_cds(self, fake_ensembl, tmp_path, caplog):
        genes = [f'GENE{n}' for n in range(sources.REST_BATCH_SIZE + 2)]
        fake_ensembl['/lookup/symbol/homo_sapiens'] = lookup_genes

        def sequences(payload):
            # The second batch (GENE50, GENE51) only holds non-coding genes, which Ensembl answers with a 404
            if 'ENST_GENE50' in payload['ids']:
                raise FileNotFoundError
            return cds_of(payload, skip={'ENST_GENE7'})
        fake_ensembl['/sequence/id?type=cds'] = sequences

        path = fetch_genes(genes, cache_dir=tmp_path)

        assert 'Skipping genes without a coding sequence: GENE50, GENE51, GENE7' in caplog.text
        with open(path) as handle:
            names = [line[1:].split('|')[0] for line in handle if line.startswith('>')]
        assert len(names) == len(genes) - 3
        assert 'GENE7' not in names

    def test_default_output_file(self, fake_ensembl):
        fake_ensembl['/lookup/symbol/mus_musculus'] = lookup_genes
        fake_ensembl['/sequence/id?type=cds'] = cds_of

        path = fetch_genes(['Trp53'], species='mouse')

        assert path == str(sources.DEFAULT_CACHE_DIR / 'ensembl_mus_musculus_selected_genes.fasta')

    def test_grch37(self, fake_ensembl, tmp_path):
        fake_ensembl['/lookup/symbol/homo_sapiens'] = lookup_genes
        fake_ensembl['/sequence/id?type=cds'] = cds_of

        path = fetch_genes(['TP53'], assembly='GRCh37', cache_dir=tmp_path)

        assert path == str(tmp_path / 'ensembl_homo_sapiens_GRCh37_selected_genes.fasta')
        assert {assembly for _, _, assembly in fake_ensembl['calls']} == {'GRCh37'}

    def test_duplicate_genes_are_retrieved_once(self, fake_ensembl, tmp_path):
        fake_ensembl['/lookup/symbol/homo_sapiens'] = lookup_genes
        fake_ensembl['/sequence/id?type=cds'] = cds_of

        fetch_genes(['TP53', 'TP53'], cache_dir=tmp_path)

        assert fake_ensembl['calls'][0][1] == {'symbols': ['TP53']}

    def test_no_genes(self):
        with pytest.raises(ValueError, match='No genes'):
            fetch_genes([])

    def test_t2t_not_available(self):
        with pytest.raises(ValueError, match='T2T'):
            fetch_genes(['TP53'], assembly='T2T')


class TestInputFile:
    def test_own_file_is_returned(self):
        assert resolve_input('my.fasta') == 'my.fasta'

    @pytest.mark.parametrize('kwargs', [{'path_to_file': 'my.fasta', 'source': 'ncbi'},
                                        {'path_to_file': 'my.fasta', 'genes': ['TP53']},
                                        {'genes': ['TP53'], 'source': 'ncbi'},
                                        {}])
    def test_invalid_combinations(self, kwargs):
        with pytest.raises(ValueError):
            resolve_input(**kwargs)

    def test_no_cache_deletes_download(self, monkeypatch):
        def fake_fetch_ccds(source, species, assembly, cache_dir):
            path = os.path.join(cache_dir, 'download.fasta')
            with open(path, mode='w') as handle:
                handle.write('>a\nATG\n')
            return path
        monkeypatch.setattr(sources, 'fetch_ccds', fake_fetch_ccds)

        with input_file(source='ncbi', cache=False) as path:
            assert os.path.exists(path)
        assert not os.path.exists(os.path.dirname(path))

    def test_genes_are_retrieved_from_ensembl(self, fake_ensembl, tmp_path):
        fake_ensembl['/lookup/symbol/homo_sapiens'] = lookup_genes
        fake_ensembl['/sequence/id?type=cds'] = cds_of

        df, _, _ = geneslider(genes=['TP53', 'MYC'], cache_dir=tmp_path, codon_of_interest='AGA',
                              initial_window_size=1, max_window_size=2)

        assert df['Gene_Symbol'].tolist() == ['TP53', 'MYC'] * 2
        assert df['Gene_ID'].tolist() == ['ENSG_TP53', 'ENSG_MYC'] * 2

    def test_genes_without_cache_are_deleted(self, fake_ensembl, monkeypatch):
        fake_ensembl['/lookup/symbol/homo_sapiens'] = lookup_genes
        fake_ensembl['/sequence/id?type=cds'] = cds_of
        # Record which file the analysis read, to check it is gone afterwards
        used = []
        real_scan = core._scan

        def record_path(path_to_file, *args):
            used.append(path_to_file)
            return real_scan(path_to_file, *args)
        monkeypatch.setattr(core, '_scan', record_path)

        df, _, _ = maxslider(genes=['TP53'], cache=False, codon_of_interest='AGA', max_window_size=1)

        assert df['Codon_Frequency'].tolist() == [1.0]
        assert not os.path.exists(used[0])

    def test_no_cache_with_cache_dir(self, tmp_path):
        with pytest.raises(ValueError, match='cache_dir'):
            with input_file(source='ncbi', cache=False, cache_dir=tmp_path):
                pass
