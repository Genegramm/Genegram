from pathlib import Path

from PIL import Image, ImageChops
from genegram import read_fasta_group, parse_rna_group

root = Path(__file__).parent.resolve()
fasta = root.parent / "data" / "seq.fasta"


def test_parse_rna_group():
    for rna_group in read_fasta_group(fasta, 50):
        for index, image_actual in parse_rna_group(rna_group):
            image_expected = Image.open(
                root.parent / "data" / "parsing" / f"{rna_group[index].description}.png"
            )

            assert ImageChops.difference(image_actual, image_expected).getbbox() is None
