from pathlib import Path

__all__ = [
    "ROOT",
    "ANALYSIS",
    "PARSING",
    "NUCLEOTIDE_TO_COLOR",
]

ROOT = Path(__file__).parent.resolve()
ANALYSIS = ROOT / "analysis"
PARSING = ROOT / "parsing"

NUCLEOTIDE_TO_COLOR = {"a": 32, "c": 64, "g": 96, "u": 128}
