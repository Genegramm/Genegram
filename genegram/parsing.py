import logging
from collections import namedtuple
from pathlib import Path
from typing import Iterator

from PIL import Image, ImageDraw
from cfpq_data import cfg_from_txt
from pyformlang.cfg import Terminal

from genegram.cfpq_pyalgo import BooleanMatrixGraph, CNF, all_pairs_reachability_matrix
from genegram.shared import ROOT

GRAMMAR = CNF.from_cfg(cfg_from_txt(ROOT / "grammar.txt"))
NUCLEOTIDE_TO_COLOR = {"a": 32, "c": 64, "g": 96, "u": 128}

RNA = namedtuple("RNA", ["description", "sequence"])

__all__ = [
    "RNA",
    "read_fasta",
    "rna_to_img",
]


def read_fasta(fasta: Path) -> Iterator[RNA]:
    logging.info(f"Read {fasta=}")
    with open(fasta, "r") as fin:
        while True:
            desc = fin.readline().strip()

            # EOF
            if not desc:
                break

            seq = fin.readline().strip()

            logging.debug(f"read_fasta():\n{desc=} \n{seq=}")

            yield RNA(desc[1:], seq.lower())


def rna_to_img(rna: RNA) -> Image:
    logging.info(f"{rna=} to image")

    bmg = BooleanMatrixGraph(matrices_size=len(rna.sequence) + 1)

    for i, nucleotide in enumerate(rna.sequence):
        bmg[Terminal(nucleotide)][i, i + 1] = True

    reachabilities = all_pairs_reachability_matrix(
        graph=bmg,
        grammar=GRAMMAR,
    )

    # create white&black 8-bit image
    img = Image.new(mode="L", size=(bmg.matrices_size - 1, bmg.matrices_size - 1))
    im_draw = ImageDraw.Draw(img)

    # draw reachabilities
    I, J, _ = reachabilities.to_lists()
    for k, i in enumerate(I):
        j = J[k]
        im_draw.line(xy=[(j - 3, i + 2), (j - 1, i)], fill=255)

    # draw letters
    for i, nucleotide in enumerate(rna.sequence):
        im_draw.point(xy=(i, i), fill=NUCLEOTIDE_TO_COLOR[nucleotide])

    logging.debug(f"rna_to_img():\n{bmg=} \n{reachabilities=} \n{img=}")

    return img
