import os
from filecmp import dircmp
from pathlib import Path

root = Path(__file__).parent.resolve()
fasta = root.parent / "data" / "seq.fasta"


def test_process(tmpdir):
    tmp = tmpdir.mkdir("tmp")
    os.system(f"python -m genegram -i {fasta} -o {tmp}")

    cmp = dircmp(tmp, root.parent / "data" / "process")

    assert cmp.left_only == []
    assert cmp.right_only == []
    assert cmp.diff_files == []
    assert cmp.funny_files == []
