from collections import namedtuple
from pathlib import Path

from pyformlang.cfg import CFG

from genegram.cfpq_pyalgo import WCNF

__all__ = [
    "ROOT",
    "GROUP_LEN",
    "NUCLEOTIDE_TO_COLOR",
    "RNA",
    "GRAMMAR",
]

GROUP_LEN = 100_000

ROOT = Path(__file__).parent.resolve()

NUCLEOTIDE_TO_COLOR = {"A": 32, "C": 64, "G": 96, "U": 128}

RNA = namedtuple("RNA", ["description", "sequence"])

GRAMMAR = WCNF(
    CFG.from_text(
        (
            """
    S -> S1
    S0 -> Any_str | Any_str S1 S0
    S2 -> a S0 u | g S0 c | u S0 a | c S0 g
    S3 -> a S2 u | g S2 c | u S2 a | c S2 g
    S4 -> a S3 u | g S3 c | u S3 a | c S3 g
    S1 -> a S1 u | u S1 a | c S1 g | g S1 c | S4
    Any -> a | u | c | g
    Any_str -> Any
    Any_str -> Any Any
    Any_str -> Any Any Any
    Any_str -> Any Any Any Any
    Any_str -> Any Any Any Any Any
    Any_str -> Any Any Any Any Any Any
    Any_str -> Any Any Any Any Any Any Any
    Any_str -> Any Any Any Any Any Any Any Any
    Any_str -> Any Any Any Any Any Any Any Any Any
    Any_str -> Any Any Any Any Any Any Any Any Any Any
    Any_str -> Any Any Any Any Any Any Any Any Any Any Any
    Any_str -> Any Any Any Any Any Any Any Any Any Any Any Any
    Any_str -> Any Any Any Any Any Any Any Any Any Any Any Any Any
    Any_str -> Any Any Any Any Any Any Any Any Any Any Any Any Any Any
    Any_str -> Any Any Any Any Any Any Any Any Any Any Any Any Any Any Any
    Any_str -> Any Any Any Any Any Any Any Any Any Any Any Any Any Any Any Any
    Any_str -> Any Any Any Any Any Any Any Any Any Any Any Any Any Any Any Any Any
    Any_str -> Any Any Any Any Any Any Any Any Any Any Any Any Any Any Any Any Any Any
    Any_str -> Any Any Any Any Any Any Any Any Any Any Any Any Any Any Any Any Any Any Any
    Any_str -> Any Any Any Any Any Any Any Any Any Any Any Any Any Any Any Any Any Any Any Any
    """
        )
    )
)
