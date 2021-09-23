import tempfile
from argparse import ArgumentParser, RawTextHelpFormatter
from pathlib import Path

from genegram import seq_fasta_to_pictures, predict, ROOT

if __name__ == "__main__":
    parser = ArgumentParser(
        description="Genegram", formatter_class=RawTextHelpFormatter
    )
    parser.add_argument(
        "--seq_fasta", required=True, type=str, help="Path to the `seq.fasta` file"
    )
    parser.add_argument(
        "--out",
        required=True,
        type=str,
        help="Path to the folder where the predictions will be saved",
    )
    parser.add_argument(
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

    seq_fasta = Path(args.seq_fasta).resolve()

    with tempfile.TemporaryDirectory() as tmp_dir:
        seq_fasta_to_pictures(seq_fasta, tmp_dir)
        predict(tmp_dir, args.out, ROOT / "weights" / f"{args.model}.h5")
