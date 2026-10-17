import logging

import pytest

from codonslider import geneslider, maxslider, slidwin


@pytest.fixture
def fasta(tmp_path):
    path = tmp_path / 'test.fasta'
    path.write_text('>g1\n'
                    'ATGAGAAGAGCC\n'
                    '>g2 gene:ENSG00000000002.3 gene_symbol:FOO\n'
                    'AGAAGAAGATTT\n')
    return str(path)


def write_fasta(tmp_path, text):
    path = tmp_path / 'sequences.fasta'
    path.write_text(text)
    return str(path)


class TestSlidwin:
    def test_finds_densest_window(self):
        result = slidwin('ATGAGAAGAGCC', 2, 1, 'AGA', 'g1')
        assert result == {'Window_Size': 2, 'Codon_Frequency': 1.0, 'Sequence': 'AGAAGA', 'Gene': 'g1'}

    def test_partial_density(self):
        result = slidwin('ATGAGAAGAGCC', 3, 1, 'AGA')
        assert result['Codon_Frequency'] == pytest.approx(2 / 3)

    def test_codon_absent(self):
        # Every window ties at 0; the last one (here the stop codon) must not be shown as a stretch
        result = slidwin('ATGCCCGGGTAA', 1, 1, 'AGA')
        assert result['Codon_Frequency'] == 0
        assert result['Sequence'] == ''

    def test_only_counts_in_frame(self):
        # 'AGA' appears out of frame (positions 1-3) but never as a codon
        result = slidwin('CAGATTCCC', 1, 1, 'AGA')
        assert result['Codon_Frequency'] == 0
        assert result['Sequence'] == ''

    def test_step_skips_windows(self):
        # With a step of 2 codons the window AGAAGA (offset 1 codon) is never visited
        assert slidwin('CCCAGAAGACCC', 2, 1, 'AGA')['Codon_Frequency'] == 1.0
        assert slidwin('CCCAGAAGACCC', 2, 2, 'AGA')['Codon_Frequency'] == 0.5

    def test_tie_keeps_last_window(self):
        result = slidwin('AGACCCCCCAGA', 2, 1, 'AGA')
        assert result['Codon_Frequency'] == 0.5
        assert result['Sequence'] == 'CCCAGA'

    def test_sequence_as_long_as_window(self):
        result = slidwin('AGAAGA', 2, 1, 'AGA')
        assert result['Codon_Frequency'] == 1.0
        assert result['Sequence'] == 'AGAAGA'

    def test_sequence_shorter_than_window(self):
        result = slidwin('AGA', 2, 1, 'AGA', 'g1')
        assert result == {'Window_Size': 2, 'Codon_Frequency': 0, 'Sequence': '', 'Gene': 'g1'}

    def test_incomplete_last_codon_is_ignored(self):
        assert slidwin('AGAAG', 1, 1, 'AGA')['Sequence'] == 'AGA'
        assert slidwin('AGAAG', 2, 1, 'AGA')['Sequence'] == ''

    def test_case_insensitive(self):
        result = slidwin('atgagaaga', 2, 1, 'aga')
        assert result['Codon_Frequency'] == 1.0
        assert result['Sequence'] == 'AGAAGA'

    def test_rna_codon(self):
        assert slidwin('ATGAGT', 1, 1, 'AGU')['Codon_Frequency'] == 1.0

    @pytest.mark.parametrize('codon', ['AG', 'AGAA', 'AGX', ''])
    def test_invalid_codon(self, codon):
        with pytest.raises(ValueError, match='not a codon'):
            slidwin('ATGAGA', 1, 1, codon)

    def test_codon_not_a_string(self):
        with pytest.raises(TypeError):
            slidwin('ATGAGA', 1, 1, None)

    @pytest.mark.parametrize('window_size, step', [(0, 1), (1, 0), (-1, 1), (1.5, 1), (True, 1)])
    def test_invalid_window_or_step(self, window_size, step):
        with pytest.raises(ValueError, match='at least 1'):
            slidwin('ATGAGA', window_size, step, 'AGA')


class TestMaxslider:
    def test_results(self, fasta):
        df, window_sizes, frequencies = maxslider(path_to_file=fasta, codon_of_interest='AGA',
                                                  initial_window_size=1, max_window_size=3)
        # No gene columns: the densest windows are often in many genes at once, and geneslider lists them
        assert list(df.columns) == ['Window_Size', 'Codon_Frequency', 'Sequence']
        assert df['Window_Size'].tolist() == [1, 2, 3]
        assert df['Codon_Frequency'].tolist() == [1.0, 1.0, 1.0]
        assert df['Sequence'].tolist() == ['AGA', 'AGAAGA', 'AGAAGAAGA']
        assert window_sizes == [1, 2, 3]
        assert frequencies == [1.0, 1.0, 1.0]

    def test_logs_highest_frequency(self, tmp_path, caplog):
        # The densest gene comes first, so the message must report its frequency, not the last gene's
        path = write_fasta(tmp_path, '>dense\nAGAAGA\n>sparse\nATGATG\n')
        with caplog.at_level(logging.INFO):
            maxslider(path_to_file=path, codon_of_interest='AGA', initial_window_size=2, max_window_size=2)
        assert 'Window size 2: highest frequency 1.000' in caplog.text
        assert 'dense' not in caplog.text.replace(path, '')

    def test_initial_window_size(self, fasta):
        df, window_sizes, _ = maxslider(path_to_file=fasta, codon_of_interest='AGA',
                                        initial_window_size=2, max_window_size=3)
        assert window_sizes == [2, 3]
        assert len(df) == 2

    def test_window_as_long_as_gene(self, tmp_path):
        path = write_fasta(tmp_path, '>a\nAGAAGA\n')
        df, _, _ = maxslider(path_to_file=path, codon_of_interest='AGA', initial_window_size=2, max_window_size=2)
        assert df.loc[0, ['Codon_Frequency', 'Sequence']].tolist() == [1.0, 'AGAAGA']

    def test_window_sizes_longer_than_every_gene_are_left_out(self, fasta):
        df, window_sizes, frequencies = maxslider(path_to_file=fasta, codon_of_interest='AGA',
                                                  initial_window_size=3, max_window_size=6)
        # Both genes are 4 codons long
        assert window_sizes == [3, 4]
        assert frequencies == [1.0, 0.75]
        assert df['Window_Size'].tolist() == [3, 4]

    def test_skips_genes_shorter_than_window(self, tmp_path):
        # The short gene has the densest stretch for a window of 2, but is too short for a window of 3
        path = write_fasta(tmp_path, '>long\nAGACCCAGATTT\n>short\nAGAAGA\n')
        df, _, _ = maxslider(path_to_file=path, codon_of_interest='AGA', initial_window_size=3, max_window_size=3)
        assert df['Codon_Frequency'].tolist() == [pytest.approx(2 / 3)]
        assert df['Sequence'].tolist() == ['AGACCCAGA']

    def test_no_gene_has_codon(self, tmp_path, caplog):
        path = write_fasta(tmp_path, '>a gene_symbol:A\nATGCCCTAA\n>b gene_symbol:B\nATGGGGTAA\n')
        with caplog.at_level(logging.INFO):
            df, _, frequencies = maxslider(path_to_file=path, codon_of_interest='AGA', max_window_size=2)
        assert frequencies == [0, 0]
        assert df['Sequence'].tolist() == ['', '']
        assert 'Window size 1: no window contains AGA' in caplog.text

    def test_empty_file(self, tmp_path):
        path = write_fasta(tmp_path, '')
        df, window_sizes, frequencies = maxslider(path_to_file=path, codon_of_interest='AGA', max_window_size=3)
        assert df.empty
        assert list(df.columns) == ['Window_Size', 'Codon_Frequency', 'Sequence']
        assert window_sizes == frequencies == []

    def test_requires_codon_and_max_window(self, fasta):
        with pytest.raises(TypeError):
            maxslider(path_to_file=fasta, max_window_size=3)
        with pytest.raises(TypeError):
            maxslider(path_to_file=fasta, codon_of_interest='AGA')

    @pytest.mark.parametrize('kwargs, message', [
        ({'codon_of_interest': 'XYZ'}, 'not a codon'),
        ({'steps': 0}, 'steps'),
        ({'initial_window_size': 0}, 'initial_window_size'),
        ({'initial_window_size': 4, 'max_window_size': 3}, 'smaller than'),
    ])
    def test_invalid_arguments(self, fasta, kwargs, message):
        arguments = {'path_to_file': fasta, 'codon_of_interest': 'AGA', 'max_window_size': 3, **kwargs}
        with pytest.raises(ValueError, match=message):
            maxslider(**arguments)

    def test_rejects_file_and_source(self, fasta):
        with pytest.raises(ValueError):
            maxslider(path_to_file=fasta, source='ncbi', codon_of_interest='AGA', max_window_size=3)


class TestGeneslider:
    def test_one_row_per_gene_and_window(self, fasta):
        df, window_sizes, frequencies = geneslider(path_to_file=fasta, codon_of_interest='AGA',
                                                   initial_window_size=1, max_window_size=3)
        assert len(df) == 6
        assert df['Window_Size'].tolist() == [1, 1, 2, 2, 3, 3]
        assert df['Gene'].tolist() == ['g1', 'g2'] * 3
        assert window_sizes == [1, 2, 3]
        assert frequencies == [1.0, 1.0, 1.0]

        largest = df[df['Window_Size'] == 3].set_index('Gene')
        assert largest.loc['g1', 'Codon_Frequency'] == pytest.approx(2 / 3)
        assert largest.loc['g2', 'Codon_Frequency'] == 1.0
        assert largest.loc['g2', 'Sequence'] == 'AGAAGAAGA'

    def test_gene_columns_from_header(self, fasta):
        df, _, _ = geneslider(path_to_file=fasta, codon_of_interest='AGA', initial_window_size=1, max_window_size=1)
        rows = df.set_index('Gene')
        assert (rows.loc['g1', 'Gene_Symbol'], rows.loc['g1', 'Gene_ID']) == ('', '')
        assert (rows.loc['g2', 'Gene_Symbol'], rows.loc['g2', 'Gene_ID']) == ('FOO', 'ENSG00000000002.3')

    def test_skips_genes_shorter_than_window(self, tmp_path):
        path = write_fasta(tmp_path, '>long\nAGAAGAAGATTT\n>short\nAGAAGA\n')
        df, _, _ = geneslider(path_to_file=path, codon_of_interest='AGA', initial_window_size=3, max_window_size=3)
        assert df['Gene'].tolist() == ['long']

    def test_window_as_long_as_gene(self, tmp_path):
        path = write_fasta(tmp_path, '>a\nAGAAGA\n')
        df, _, _ = geneslider(path_to_file=path, codon_of_interest='AGA', initial_window_size=2, max_window_size=2)
        assert df.loc[0, ['Codon_Frequency', 'Gene']].tolist() == [1.0, 'a']

    def test_gene_without_codon(self, tmp_path):
        # The gene keeps its row and name, but has no window to show
        path = write_fasta(tmp_path, '>with gene_symbol:W\nAGAATGTAA\n>without gene_symbol:WO\nATGCCCTAA\n')
        df, _, frequencies = geneslider(path_to_file=path, codon_of_interest='AGA', max_window_size=1)
        rows = df.set_index('Gene')
        assert rows.loc['with', ['Codon_Frequency', 'Sequence']].tolist() == [1.0, 'AGA']
        assert rows.loc['without', ['Codon_Frequency', 'Sequence', 'Gene_Symbol']].tolist() == [0, '', 'WO']
        assert frequencies == [1.0]

    def test_window_longer_than_every_gene(self, fasta, caplog):
        with caplog.at_level(logging.INFO):
            df, window_sizes, frequencies = geneslider(path_to_file=fasta, codon_of_interest='AGA',
                                                       initial_window_size=5, max_window_size=5)
        assert df.empty
        assert window_sizes == []
        assert frequencies == []
        assert 'Window size 5: no sequence is that long' in caplog.text

    def test_requires_codon_and_max_window(self, fasta):
        with pytest.raises(TypeError):
            geneslider(path_to_file=fasta, max_window_size=3)
        with pytest.raises(TypeError):
            geneslider(path_to_file=fasta, codon_of_interest='AGA')
