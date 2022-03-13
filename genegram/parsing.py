"""FASTA parsing module"""
from pathlib import Path
from typing import List, Tuple

import numpy as np

from genegram.cfpq_pyalgo import BooleanMatrixGraph, all_pairs_reachability_matrix
from genegram.shared import GROUP_LEN, GRAMMAR, RNA, NUCLEOTIDE_TO_COLOR

__all__ = [
    "read_fasta_single",
    "read_fasta_group",
    "parse_rna_single",
    "parse_rna_group",
]


def read_fasta_single(fasta: Path) -> List[RNA]:
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
    rna: List[RNA]
        RNA sequences from a FASTA file
    """
    with open(fasta, "r") as fin:
        data = fin.readlines()
    return [
        RNA(data[i].strip()[1:], data[i + 1].strip()) for i in range(0, len(data), 2)
    ]


def read_fasta_group(fasta: Path, limit: int = GROUP_LEN) -> List[List[RNA]]:
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
    rna: List[List[RNA]]
        Grouped RNA sequences from a FASTA file
    """
    with open(fasta, "r") as fin:
        data = fin.readlines()
    result = []
    cur_group = []
    cur_len = 0

    for i in range(0, len(data), 2):
        rna = RNA(data[i].strip()[1:], data[i + 1].strip())
        cur_group.append(rna)
        cur_len += 1 + len(rna.sequence)

        if abs(limit - cur_len) < abs(limit - (cur_len + 1 + len(rna.sequence))):
            result.append(cur_group)
            cur_group = []
            cur_len = 0

    if cur_group:
        result.append(cur_group)

    return result


def parse_rna_single(rna: RNA) -> np.ndarray:
    """Parse RNA `rna` with grammar `GRAMMAR`

    Parameters
    ----------
    rna: RNA

    Returns
    -------
    image: np.ndarray
        The result of parsing presented as an image
    """
    n = len(rna.sequence)
    bmg = BooleanMatrixGraph(n + 1)
    # in grammar, nucleotides must be lowercase
    for i, nucleotide in enumerate(rna.sequence.lower()):
        bmg[nucleotide][i, i + 1] = True

    reachabilities = all_pairs_reachability_matrix(bmg, GRAMMAR)

    image = np.zeros((n, n), dtype=np.uint8)

    # draw pairings
    I, J, _ = reachabilities.to_lists()
    for i, j in zip(I, J):
        image[i + 2, j - 3] = 255
        image[i + 1, j - 2] = 255
        image[i, j - 1] = 255

    # draw nucleotides
    for i, nucleotide in enumerate(rna.sequence):
        image[i, i] = NUCLEOTIDE_TO_COLOR[nucleotide]

    return image


def parse_rna_group(rna_group: List[RNA]) -> List[Tuple[int, np.ndarray]]:
    """Parse a group of RNA sequences `rna_group` with grammar `GRAMMAR`

    Parameters
    ----------
    rna_group: List[RNA]
        A group of RNA sequences

    Returns
    -------
    index, image: Tuple[int, np.ndarray]
        Iterator over the result of parsing
    """
    result = []
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
        image = np.zeros((n, n), dtype=np.uint8)

        # draw pairings
        I, J, _ = reachabilities[
            prefix : (prefix + n), prefix : (prefix + n)
        ].to_lists()
        for i, j in zip(I, J):
            image[i + 2, j - 3] = 255
            image[i + 1, j - 2] = 255
            image[i, j - 1] = 255

        # draw nucleotides
        for i, nucleotide in enumerate(seq):
            image[i, i] = NUCLEOTIDE_TO_COLOR[nucleotide]

        prefix += n + 1
        result.append((index, image))

    return result
