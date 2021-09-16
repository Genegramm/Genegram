import tempfile
from argparse import ArgumentParser
from pathlib import Path

from genegram import seq_fasta_to_pictures, predict

if __name__ == "__main__":
    parser = ArgumentParser(description="Genegram")
    parser.add_argument("--seq_fasta", required=True, type=str, help="input data path")
    parser.add_argument(
        "--predict",
        required=True,
        type=str,
        help="output data path",
    )

    args = parser.parse_args()

    seq_fasta = Path(args.seq_fasta).resolve()

    with tempfile.TemporaryDirectory() as tmp_dir:
        seq_fasta_to_pictures(seq_fasta, tmp_dir)
        predict(tmp_dir, args.predict)
