from pathlib import Path

import numpy as np
from PIL import Image

from genegram import read_fasta_single, parse_rna_single

root = Path(__file__).parent.resolve()
fasta = root.parent / "data" / "seq.fasta"


def test_parse_rna_single():
    for rna in read_fasta_single(fasta):
        image_actual = parse_rna_single(rna)
        image_expected = np.array(
            Image.open(root.parent / "data" / "parsing" / f"{rna.description}.png")
        )

        assert np.array_equal(image_actual, image_expected)
