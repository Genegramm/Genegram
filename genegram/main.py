"""FASTA processing module"""
from pathlib import Path

from genegram.parsing import (
    read_fasta_group,
    parse_rna_group,
    read_fasta_single,
    parse_rna_single,
)
from genegram.predict import setup_model, clear_session, predict
from genegram.shared import ROOT
from genegram.utils import binarize_image, create_connectivity_table, remove_multiplets

__all__ = [
    "process_fasta_single",
    "process_fasta_group",
]


def process_fasta_single(
    fasta: Path,
    out: Path,
    weights: str = "main",
    bin_coeff: float = 0.6,
) -> Path:
    """Process Secondary Structure for each RNA from `fasta` file using single method

    Parameters
    ----------
    fasta: Path
        The path to FASTA input file
    out:
        The path where the resulting connectivity tables will be saved
    weights:
        Type of predictive model weights
    bin_coeff:
        Binarization coefficient

    Returns
    -------
    out: Path
        The path where the resulting connectivity tables will be saved
    """
    if not out.exists():
        out.mkdir(parents=True, exist_ok=True)

    model = setup_model(ROOT / "weights" / f"{weights}.h5")

    for rna in read_fasta_single(fasta):
        image = parse_rna_single(rna)
        prediction = predict(image, model)
        prediction_cleaned = remove_multiplets(prediction)
        pred_bin = binarize_image(prediction_cleaned, bin_coeff)
        ct = create_connectivity_table(pred_bin, *rna)

        target_path = out / f"{rna.description}.ct"
        with open(target_path, "w") as fout:
            fout.write(ct)

    clear_session()

    return Path(out).resolve()


def process_fasta_group(
    fasta: Path,
    out: Path,
    weights: str = "main",
    bin_coeff: float = 0.6,
) -> Path:
    """Process Secondary Structure for each RNA from `fasta` file using group method

    Parameters
    ----------
    fasta: Path
        The path to FASTA input file
    out:
        The path where the resulting connectivity tables will be saved
    weights:
        Type of predictive model weights
    bin_coeff:
        Binarization coefficient

    Returns
    -------
    out: Path
        The path where the resulting connectivity tables will be saved
    """
    if not out.exists():
        out.mkdir(parents=True, exist_ok=True)

    model = setup_model(ROOT / "weights" / f"{weights}.h5")

    for rna_group in read_fasta_group(fasta):
        for index, image in parse_rna_group(rna_group):
            prediction = predict(image, model)
            prediction_cleaned = remove_multiplets(prediction)
            pred_bin = binarize_image(prediction_cleaned, bin_coeff)
            ct = create_connectivity_table(pred_bin, *rna_group[index])

            target_path = out / f"{rna_group[index].description}.ct"
            with open(target_path, "w") as fout:
                fout.write(ct)

    clear_session()

    return Path(out).resolve()
