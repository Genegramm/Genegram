from pathlib import Path

from PIL import Image, ImageChops
from genegram import read_fasta_single, parse_rna_single

root = Path(__file__).parent.resolve()
fasta = root.parent / "data" / "seq.fasta"


def test_parse_rna_single():
    for rna in read_fasta_single(fasta):
        image_actual = parse_rna_single(rna)
        image_expected = Image.open(
            root.parent / "data" / "parsing" / f"{rna.description}.png"
        )

        assert ImageChops.difference(image_actual, image_expected).getbbox() is None
