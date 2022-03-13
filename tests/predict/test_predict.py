from glob import glob
from pathlib import Path

import numpy as np
from PIL import Image, ImageChops

from genegram import setup_model, predict, clear_session, __path__ as genegram_path

root = Path(__file__).parent.resolve()


def test_predict():
    model = setup_model(Path(genegram_path[0]) / "weights" / "main.h5")
    for image_path in glob(str(root.parent / "data" / "parsing" / "*")):
        image = Image.open(image_path)
        actual_prediction = Image.fromarray(predict(np.array(image), model))
        expected_prediction = Image.open(
            root.parent / "data" / "predict" / f"{Path(image_path).stem}.png"
        )

        assert (
            ImageChops.difference(actual_prediction, expected_prediction).getbbox()
            is None
        )
    clear_session()
