from argparse import ArgumentParser, RawTextHelpFormatter
from pathlib import Path

from genegram.main import process_fasta

if __name__ == "__main__":
    parser = ArgumentParser(
        description="Genegram", formatter_class=RawTextHelpFormatter
    )
    parser.add_argument(
        "-i", "--inp", required=True, type=str, help="Path to the FASTA file"
    )
    parser.add_argument(
        "-o",
        "--out",
        required=True,
        type=str,
        help="Path to the folder where the predictions will be saved",
    )
    parser.add_argument(
        "-m",
        "--model",
        required=False,
        type=str,
        choices=["main", "mps", "pks"],
        default="main",
        help=(
            "Type of the model to be used:"
            "\nmain -- The default model, the best on average"
            "\nmps -- Multiplet prediction model"
            "\npks -- Pseudoknots prediction model"
        ),
    )

    args = parser.parse_args()

    process_fasta(
        fasta=Path(args.inp).resolve(),
        out=Path(args.out).resolve(),
        weights=args.model,
    )
