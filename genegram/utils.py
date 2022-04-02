"""Post-processing utilities"""
__all__ = [
    "binarize_image",
    "create_connectivity_table",
    "remove_multiplets",
]

import numpy as np


def create_connectivity_table(image: np.ndarray, meta: str, seq: str) -> str:
    """Creates a Connectivity Table from an RNA Secondary Structure
    represented as an image

    Parameters
    ----------
    image: np.ndarray
        RNA Secondary Structure
    meta: str
        RNA description
    seq: str
        RNA sequence

    Returns
    -------
    ct: str
        Connectivity Table
    """
    size = len(image)
    # ct = "  " + str(size) + " " + meta + "\n"
    ct = f"  {size} {meta}\n"
    for i in range(size):
        pair = 0
        for j in range(size):
            if image[i][j] == 255 or image[j][i] == 255:
                pair = j + 1
        ct += f"    {i + 1} {seq[i]}       {i}    {i + 2}  {pair}    {i + 1}\n"
        # ct += (
        #     "    "
        #     + str(i + 1)
        #     + " "
        #     + seq[i]
        #     + "       "
        #     + str(i)
        #     + "    "
        #     + str(i + 2)
        #     + "  "
        #     + str(pair)
        #     + "    "
        #     + str(i + 1)
        #     + "\n"
        # )
    return ct


def remove_multiplets(image: np.ndarray) -> np.ndarray:
    """Remove multiplets from network prediction of RNA Secondary Structure

    Parameters
    ----------
    image: np.ndarray
        RNA Secondary Structure

    Returns
    -------
    image: np.ndarray
        RNA Secondary Structure with removed multiplets
    """

    def get_multiplets(i0, j0, image):
        mps = []
        size = len(image)
        for i in range(size):
            if image[i0, i] == 255 and (i0, i) != (i0, j0) and i0 <= i:
                mps.append((i0, i))
            if image[j0, i] == 255 and (j0, i) != (i0, j0) and j0 <= i:
                mps.append((j0, i))
            if image[i, i0] == 255 and (i, i0) != (i0, j0) and i <= i0:
                mps.append((i, i0))
            if image[i, j0] == 255 and (i, j0) != (i0, j0) and i <= j0:
                mps.append((i, j0))
        return list(set(mps))

    def get_stem_len(i0, j0, image):
        size = len(image)
        stem_len = 1
        i, j = i0 + 1, j0 - 1
        while (
            i < len(image)
            and j >= 0
            and image[i][j] == 255
            and len(get_multiplets(i, j, image)) == 0
        ):
            stem_len += 1
            i += 1
            j -= 1
        i, j = i0 - 1, j0 + 1
        while (
            i >= 0
            and j < len(image)
            and image[i][j] == 255
            and len(get_multiplets(i, j, image)) == 0
        ):
            stem_len += 1
            i -= 1
            j += 1
        return stem_len

    size = len(image)
    mps_nums = dict()
    for i in range(size):
        for j in range(i + 1, size):
            if image[i][j] == 255:
                mps = get_multiplets(i, j, image)
                if len(mps) > 0:
                    mps_nums[(i, j)] = len(mps)
    mps_nums = {k: v for k, v in sorted(mps_nums.items(), key=lambda item: -item[1])}
    to_delete = []
    while len(mps_nums) > 0:
        for el in to_delete:
            del mps_nums[el]
        to_delete = []
        for (i, j) in mps_nums.keys():
            if not (i, j) in to_delete:
                mps = get_multiplets(i, j, image)
                if len(mps) > 0:
                    min_stem_len = get_stem_len(i, j, image)
                    i_pick, j_pick = i, j
                    for (i0, j0) in mps:
                        stem_len = get_stem_len(i0, j0, image)
                        if stem_len < min_stem_len:
                            i_pick, j_pick = i0, j0
                            min_stem_len = stem_len
                    image[i_pick][j_pick] = 0
                    to_delete.append((i_pick, j_pick))
                else:
                    to_delete.append((i, j))
    return image


def binarize_image(image: np.ndarray, bin_coeff: float = 0.6) -> np.ndarray:
    """Set each gray pixel of network prediction `image` to black/white
    according to threshold coeff

    Parameters
    ----------
    image: np.ndarray
        Network prediction image
    bin_coeff: float
        Binarization coefficient

    Returns
    -------
    im: np.ndarray
        Binarized image
    """
    im = image.copy()
    size = len(image)
    for i in range(size):
        for j in range(size):
            if i != j:
                if im[i][j] > 255 * bin_coeff:
                    im[i][j] = 255
                else:
                    im[i][j] = 0
    return im
