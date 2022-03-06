from pathlib import Path
from typing import Iterator, List

from PIL import Image, ImageDraw
from pyformlang.cfg import Terminal

from genegram.cfpq_pyalgo import BooleanMatrixGraph, all_pairs_reachability_matrix
from genegram.shared import GROUP_LEN, GRAMMAR, RNA, NUCLEOTIDE_TO_COLOR

__all__ = [
    "read_fasta",
    "create_images",
]


def read_fasta(fasta: Path, limit: int = GROUP_LEN) -> Iterator[List[RNA]]:
    cur_group = []
    cur_len = 0

    with open(fasta, "r") as fin:
        while True:
            desc = fin.readline().strip()[1:]

            # EOF reached
            if not desc:
                break

            seq = fin.readline().strip()

            rna = RNA(desc, seq)

            if abs(limit - cur_len) > abs(limit - (cur_len + 1 + len(seq))):
                cur_group.append(rna)
                cur_len += 1 + len(seq)
            else:
                yield cur_group
                cur_group = [rna]
                cur_len = 1 + len(seq)

    if cur_group:
        yield cur_group


def create_images(rna_group: List[RNA]):
    # in grammar, terminals must be lowercase
    glued_rna = "$".join((seq.lower() for _, seq in rna_group))
    bmg = BooleanMatrixGraph(matrices_size=len(glued_rna) + 1)

    for i, nucleotide in enumerate(glued_rna):
        bmg[Terminal(nucleotide)][i, i + 1] = True

    reachabilities = all_pairs_reachability_matrix(graph=bmg, grammar=GRAMMAR)

    prefix = 0
    for index, (desc, seq) in enumerate(rna_group):
        n = len(seq)

        # create white&black 8-bit image
        image = Image.new(mode="L", size=(n, n))
        image_draw = ImageDraw.Draw(image)

        # draw pairings
        I, J, _ = reachabilities[
            prefix : (prefix + n), prefix : (prefix + n)
        ].to_lists()
        for i, j in zip(I, J):
            image_draw.line(xy=[(j - 3, i + 2), (j - 1, i)], fill=255)

        # draw nucleotides
        for i, nucleotide in enumerate(seq):
            image_draw.point(xy=(i, i), fill=NUCLEOTIDE_TO_COLOR[nucleotide])

        prefix += n + 1
        yield index, image
