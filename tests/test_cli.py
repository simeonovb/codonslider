import io
import runpy
import sys

import pandas as pd
import pytest

from codonslider import __version__
from codonslider.cli import main


@pytest.fixture
def fasta(tmp_path):
    path = tmp_path / 'test.fasta'
    path.write_text('>g1\nATGAGAAGAGCC\n>g2\nAGAAGAAGATTT\n')
    return str(path)


def test_max_writes_csv(fasta, tmp_path):
    output = tmp_path / 'out.csv'
    assert main(['max', '--codon', 'aga', '--max-window', '3', '--fasta', fasta, '-o', str(output)]) == 0
    df = pd.read_csv(output)
    assert df['Window_Size'].tolist() == [1, 2, 3]
    assert df['Sequence'].tolist() == ['AGA', 'AGAAGA', 'AGAAGAAGA']


def test_gene_writes_to_stdout(fasta, capsys):
    assert main(['gene', '-c', 'AGA', '--min-window', '2', '--max-window', '3', '--step', '1', '-f', fasta, '-q']) == 0
    df = pd.read_csv(io.StringIO(capsys.readouterr().out))
    assert len(df) == 4
    assert df['Gene'].tolist() == ['g1', 'g2'] * 2


def test_invalid_codon(fasta, capsys):
    assert main(['max', '--codon', 'XYZ', '--max-window', '3', '--fasta', fasta]) == 1
    assert "error: 'XYZ' is not a codon" in capsys.readouterr().err


def test_missing_file(tmp_path, capsys):
    assert main(['max', '--codon', 'AGA', '--max-window', '3', '--fasta', str(tmp_path / 'missing.fasta')]) == 1
    assert 'missing.fasta' in capsys.readouterr().err


def test_sequences_required(capsys):
    with pytest.raises(SystemExit) as exit_info:
        main(['max', '--codon', 'AGA', '--max-window', '3'])
    assert exit_info.value.code == 2


def test_only_one_sequence_source(fasta):
    with pytest.raises(SystemExit):
        main(['max', '--codon', 'AGA', '--max-window', '3', '--fasta', fasta, '--source', 'ncbi'])


def test_version(capsys):
    with pytest.raises(SystemExit):
        main(['--version'])
    assert capsys.readouterr().out.strip() == f'codonslider {__version__}'


def test_python_m(monkeypatch, capsys):
    monkeypatch.setattr(sys, 'argv', ['codonslider', '--version'])
    with pytest.raises(SystemExit):
        runpy.run_module('codonslider', run_name='__main__')
    assert __version__ in capsys.readouterr().out
