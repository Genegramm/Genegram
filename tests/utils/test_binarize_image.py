from glob import glob
from pathlib import Path

from PIL import Image, ImageChops

from genegram import (
    setup_model,
    predict,
    binarize_image,
    remove_multiplets,
    clear_session,
    __path__ as genegram_path,
)

root = Path(__file__).parent.resolve()


def test_binarize_image():
    model = setup_model(Path(genegram_path[0]) / "weights" / "main.h5")
    for image_path in glob(str(root.parent / "data" / "parsing" / "*")):
        image = Image.open(image_path)
        prediction = predict(image, model)
        pred_bin = binarize_image(prediction)

        actual_binarized = Image.fromarray(pred_bin)
        expected_binarized = Image.open(
            root.parent / "data" / "binarize" / f"{Path(image_path).stem}.png"
        )

        assert (
            ImageChops.difference(actual_binarized, expected_binarized).getbbox()
            is None
        )
    clear_session()
