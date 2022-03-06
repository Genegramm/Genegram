from pathlib import Path

from genegram.parsing import read_fasta, create_images
from genegram.predict import setup_model, clear_session, predict
from genegram.shared import ROOT
from genegram.utils import binarize_image, create_connectivity_table

__all__ = [
    "process_fasta",
]


def process_fasta(
    fasta: Path,
    out: Path,
    weights: str = "main",
    bin_coeff: float = 0.6,
):
    if not out.exists():
        out.mkdir(parents=True, exist_ok=True)

    model = setup_model(ROOT / "weights" / f"{weights}.h5")

    for rna_group in read_fasta(fasta):
        for index, image in create_images(rna_group):
            prediction = predict(image, model)
            pred_bin = binarize_image(prediction, bin_coeff)
            ct = create_connectivity_table(pred_bin, *rna_group[index])

            target_path = out / f"{rna_group[index].description}.ct"
            with open(target_path, "w") as fout:
                fout.write(ct)

    clear_session()
