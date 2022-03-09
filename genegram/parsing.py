"""FASTA parsing module"""
from pathlib import Path
from typing import Iterator, List, Tuple

from PIL import Image, ImageDraw

from genegram.cfpq_pyalgo import BooleanMatrixGraph, all_pairs_reachability_matrix
from genegram.shared import GROUP_LEN, GRAMMAR, RNA, NUCLEOTIDE_TO_COLOR

__all__ = [
    "read_fasta_single",
    "read_fasta_group",
    "parse_rna_single",
    "parse_rna_group",
]


def read_fasta_single(fasta: Path) -> Iterator[RNA]:
    """Read FASTA file from `fasta` path.
    The FASTA file must be in the following format:
    >`RNA description`
    `RNA sequence`
    >`RNA description`
    `RNA sequence`
    ...
    >`RNA description`
    `RNA sequence`

    Parameters
    ----------
    fasta: Path
        Path to the FASTA file

    Returns
    -------
    rna: RNA
        Iterator over RNA sequences from a file
    """
    with open(fasta, "r") as fin:
        while True:
            desc = fin.readline().strip()[1:]

            # EOF reached
            if not desc:
                break

            seq = fin.readline().strip()

            yield RNA(desc, seq)


def read_fasta_group(fasta: Path, limit: int = GROUP_LEN) -> Iterator[List[RNA]]:
    """Read FASTA file from `fasta` path.
    The FASTA file must be in the following format:
    >`RNA description`
    `RNA sequence`
    >`RNA description`
    `RNA sequence`
    ...
    >`RNA description`
    `RNA sequence`

    Parameters
    ----------
    fasta: Path
        Path to the FASTA file
    limit: int
        Limit on the total size of a group of RNA sequences

    Returns
    -------
    rna_group:
        Iterator over group of RNA sequences from a file
    """
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
                if cur_group:
                    yield cur_group
                cur_group = [rna]
                cur_len = 1 + len(seq)

    if cur_group:
        yield cur_group


def parse_rna_single(rna: RNA) -> Image.Image:
    """Parse RNA `rna` with grammar `GRAMMAR`

    Parameters
    ----------
    rna: RNA

    Returns
    -------
    image: PIL.Image.Image
        The result of parsing presented as an image
    """
    n = len(rna.sequence)
    bmg = BooleanMatrixGraph(n + 1)
    # in grammar, nucleotides must be lowercase
    for i, nucleotide in enumerate(rna.sequence.lower()):
        bmg[nucleotide][i, i + 1] = True

    reachabilities = all_pairs_reachability_matrix(bmg, GRAMMAR)

    image = Image.new(mode="L", size=(n, n))
    image_draw = ImageDraw.Draw(image)

    # draw pairings
    I, J, _ = reachabilities.to_lists()
    for i, j in zip(I, J):
        image_draw.line(xy=[(j - 3, i + 2), (j - 1, i)], fill=255)

    # draw nucleotides
    for i, nucleotide in enumerate(rna.sequence):
        image_draw.point(xy=(i, i), fill=NUCLEOTIDE_TO_COLOR[nucleotide])

    return image


def parse_rna_group(rna_group: List[RNA]) -> Iterator[Tuple[int, Image.Image]]:
    """Parse a group of RNA sequences `rna_group` with grammar `GRAMMAR`

    Parameters
    ----------
    rna_group: List[RNA]
        A group of RNA sequences

    Returns
    -------
    index, image: Tuple[int, PIL.Image.Image]
        Iterator over the result of parsing
    """
    # in grammar, nucleotides must be lowercase
    glued_rna = "$".join((seq.lower() for _, seq in rna_group))
    bmg = BooleanMatrixGraph(len(glued_rna) + 1)
    for i, nucleotide in enumerate(glued_rna):
        bmg[nucleotide][i, i + 1] = True

    reachabilities = all_pairs_reachability_matrix(bmg, GRAMMAR)

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
