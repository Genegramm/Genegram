import logging
from argparse import ArgumentParser, RawTextHelpFormatter
from pathlib import Path

from genegram.parsing import read_fasta
from genegram.predict import rna_predict, setup_model, clear_session
from genegram.shared import ROOT

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
    parser.add_argument(
        "-l",
        "--log",
        required=False,
        type=str,
        choices=["INFO", "WARNING", "ERROR", "CRITICAL", "DEBUG"],
        default="INFO",
        help=(
            "Type of the logging level to be used:"
            "\nINFO -- Confirmation that things are working as expected"
            "\nWARNING -- An indication that something unexpected happened, the software is still working as expected"
            "\nERROR -- Due to a more serious problem, the software has not been able to perform some function"
            "\nCRITICAL -- A serious error, indicating that the program itself may be unable to continue running"
            "\nDEBUG -- Detailed information, typically of interest only when diagnosing problems"
        ),
    )

    args = parser.parse_args()

    logging.basicConfig(
        level=args.log,
        format="[%(asctime)s]>%(levelname)s>%(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    logging.info(f"Parse {args=}")

    out = Path(args.out).resolve()
    if not out.exists():
        out.mkdir(parents=True, exist_ok=True)
        logging.info(f"Create {out} dir")

    model = setup_model(ROOT / "weights" / f"{args.model}.h5")

    for rna in read_fasta(Path(args.inp).resolve()):
        pred = rna_predict(rna, model)

        target_path = out / f"{rna.description}.ct"

        with open(target_path, "w") as fout:
            fout.write(pred.ct)
            logging.info(
                f"Save {rna=} secondary structure connectivity table to {target_path=}"
            )

    clear_session()
