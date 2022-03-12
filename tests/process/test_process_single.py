from filecmp import dircmp
from pathlib import Path

from genegram import process_fasta_single

root = Path(__file__).parent.resolve()
fasta = root.parent / "data" / "seq.fasta"


def test_process_fasta_single(tmpdir):
    tmp = tmpdir.mkdir("tmp")
    process_fasta_single(fasta, tmp)

    cmp = dircmp(tmp, root.parent / "data" / "process")

    assert cmp.left_only == []
    assert cmp.right_only == []
    assert cmp.diff_files == []
    assert cmp.funny_files == []
