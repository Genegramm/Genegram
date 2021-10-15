import logging
from collections import namedtuple
from pathlib import Path

import numpy as np
import tensorflow as tf
from keras import backend as K
from keras import regularizers
from keras.layers import Dropout, Conv2D, Input, Activation, BatchNormalization, add
from keras.models import Model
from tensorflow.keras.layers import Layer

from genegram.parsing import RNA, rna_to_img

__all__ = [
    "PredictionContext",
    "rna_predict",
    "img_predict",
    "img_binarize",
    "img_to_ct",
    "remove_multiplets",
    "setup_model",
    "clear_session",
]

PredictionContext = namedtuple(
    "PredictionContext",
    [
        "rna",
        "img_rna",
        "img_pred",
        "ct",
    ],
)


def clear_session():
    logging.info("Clear TensorFlow session")
    K.clear_session()


def img_to_ct(img, meta, seq):
    logging.info(f"Create Connectivity Table for RNA({meta}, {seq})")

    size = len(img)
    ct = "  " + str(size) + " " + meta + "\n"
    for i in range(size):
        pair = 0
        for j in range(size):
            if img[i][j] == 255 or img[j][i] == 255:
                pair = j + 1
        ct += (
            "    "
            + str(i + 1)
            + " "
            + seq[i]
            + "       "
            + str(i)
            + "    "
            + str(i + 2)
            + "  "
            + str(pair)
            + "    "
            + str(i + 1)
            + "\n"
        )

    logging.debug(f"img_to_ct():\n{img=} \n{meta=} \n{seq=} \n{ct=}")

    return ct


def remove_multiplets(img):
    def get_multiplets(i0, j0, img):
        mps = []
        size = len(img)
        for i in range(size):
            if img[i0, i] == 255 and (i0, i) != (i0, j0) and i0 <= i:
                mps.append((i0, i))
            if img[j0, i] == 255 and (j0, i) != (i0, j0) and j0 <= i:
                mps.append((j0, i))
            if img[i, i0] == 255 and (i, i0) != (i0, j0) and i <= i0:
                mps.append((i, i0))
            if img[i, j0] == 255 and (i, j0) != (i0, j0) and i <= j0:
                mps.append((i, j0))
        return list(set(mps))

    def get_stem_len(i0, j0, img):
        size = len(img)
        stem_len = 1
        i, j = i0 + 1, j0 - 1
        while (
            i < len(img)
            and j >= 0
            and img[i][j] == 255
            and len(get_multiplets(i, j, img)) == 0
        ):
            stem_len += 1
            i += 1
            j -= 1
        i, j = i0 - 1, j0 + 1
        while (
            i >= 0
            and j < len(img)
            and img[i][j] == 255
            and len(get_multiplets(i, j, img)) == 0
        ):
            stem_len += 1
            i -= 1
            j += 1
        return stem_len

    size = len(img)
    mps_nums = dict()
    for i in range(size):
        for j in range(i + 1, size):
            if img[i][j] == 255:
                mps = get_multiplets(i, j, img)
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
                mps = get_multiplets(i, j, img)
                if len(mps) > 0:
                    min_stem_len = get_stem_len(i, j, img)
                    i_pick, j_pick = i, j
                    for (i0, j0) in mps:
                        stem_len = get_stem_len(i0, j0, img)
                        if stem_len < min_stem_len:
                            i_pick, j_pick = i0, j0
                            min_stem_len = stem_len
                    img[i_pick][j_pick] = 0
                    to_delete.append((i_pick, j_pick))
                else:
                    to_delete.append((i, j))
    return img


# set each gray pixel of network prediction to black/white
# according to threshold coeff
def img_binarize(img, coeff=0.6):
    logging.info(f"Binarize image")

    im = img.copy()
    size = len(img)
    for i in range(size):
        for j in range(size):
            if i != j:
                if im[i][j] > 255 * coeff:
                    im[i][j] = 255
                else:
                    im[i][j] = 0

    logging.debug(f"img_binarize():\n{img=} \n{coeff=} \n{im=}")

    return im


# Model definition functions

# layer that for inputs i1, i2, i3, i4
# returns (w1*i1 + w2*i2 + w3*i3 + w4*i4) / 4, where wi are trainable coeffs
class WeightedSum(Layer):
    def __init__(self):
        super(WeightedSum, self).__init__()
        w_init = tf.keras.initializers.Ones()
        self.w1 = tf.Variable(
            initial_value=w_init(shape=(), dtype="float32"), trainable=False
        )
        self.w2 = tf.Variable(
            initial_value=w_init(shape=(), dtype="float32"), trainable=False
        )
        self.w3 = tf.Variable(
            initial_value=w_init(shape=(), dtype="float32"), trainable=False
        )
        self.w4 = tf.Variable(
            initial_value=w_init(shape=(), dtype="float32"), trainable=False
        )

    def call(self, input):
        input1, input2, input3, input4 = input
        return (
            tf.multiply(input1, self.w1)
            + tf.multiply(input2, self.w2)
            + tf.multiply(input3, self.w3)
            + tf.multiply(input4, self.w4)
        ) / 4


# residual unit definition, classical structure in ML
def res_unit(inputs, filters, kernels, activ, bn=False, dr=0):
    x = inputs
    if bn:
        x = BatchNormalization()(x)
    for i in range(len(filters)):
        x = Activation(activ)(x)
        x = Conv2D(
            filters[i],
            kernel_size=kernels[i],
            activity_regularizer=regularizers.l2(1e-10),
            padding="same",
        )(x)
        if dr > 0:
            x = Dropout(0.1)(x)
    return x


# residual network definition: repeating res units
# + skip sonnections(add layers)
def res_network(inputs, units_num, filters, kernels, activ="relu", bn=False, dr=0):
    x = res_unit(inputs, filters, kernels, activ, bn, dr)
    x = add([x, inputs])
    x = Activation(activ)(x)
    for i in range(units_num - 1):
        y = res_unit(x, filters, kernels, activ, bn, dr)
        x = add([x, y])
        x = Activation(activ)(x)
    outputs = x

    model = Model(inputs=inputs, outputs=outputs)
    return model


# model that combines several residual networks
def parallel_res_network(
    blocks_num, units_num, filters, kernels, activ="relu", bn=False, dr=0
):
    inputs = Input(shape=(None, None, 1))
    all_outputs = []
    # construct several resnets of same shape
    for i in range(blocks_num):
        model = res_network(inputs, units_num, filters, kernels, activ, bn, dr)
        all_outputs.append(model.output)

    # combine their outputs by weighted sum
    x = WeightedSum()(all_outputs)

    # add final res unit to the end of network
    y = res_unit(x, filters, kernels, activ, bn, dr)
    x = add([x, y])
    x = Activation(activ)(x)
    outputs = x
    model = Model(inputs=inputs, outputs=outputs)
    return model


def setup_model(weights: Path):
    # setup session
    tf_session_config = tf.compat.v1.ConfigProto(
        intra_op_parallelism_threads=1, inter_op_parallelism_threads=1
    )
    tf_session = tf.compat.v1.Session(
        graph=tf.compat.v1.get_default_graph(), config=tf_session_config
    )
    K.set_session(tf_session)

    logging.info("Setup TensorFlow session")

    # setup model
    model = parallel_res_network(
        blocks_num=4,
        units_num=5,
        filters=[12, 10, 8, 6, 1],
        kernels=[13, 11, 9, 7, 5],
        activ="relu",
        bn=False,
        dr=0.1,
    )
    model.load_weights(weights)

    logging.info(f"Setup model")
    logging.debug(f"setup_model():\n{weights=} \n{model=}")

    return model


# get model output on input image pil_img
def img_predict(pil_img, model):
    logging.info("Predict image")

    img = np.array(pil_img)
    n = len(img)
    inp = np.array([img]).reshape((1, n, n, 1))

    pred = model.predict(inp)

    res = np.full(shape=(n, n), fill_value=255, dtype=np.uint8)
    for i in range(n):
        for j in range(n):
            res[i][j] = min(pred[0][i][j][0], 255)

    logging.debug(
        f"img_predict():\n{pil_img=} \n{model=} \n{img=} \n{n=} \n{inp=} \n{pred=} \n{res=}"
    )

    return res


def rna_predict(rna: RNA, model):
    logging.info("Predict RNA")

    img_rna = rna_to_img(rna)

    pred = img_predict(img_rna, model)

    pred_bin = img_binarize(pred)

    ct = img_to_ct(
        pred_bin,
        rna.description,
        rna.sequence.upper(),
    )

    logging.debug(
        f"rna_predict():\n{rna=} \n{model=} \n{img_rna=} \n{pred=} \n{pred_bin=} \n{ct=}"
    )

    return PredictionContext(
        rna=rna,
        img_rna=img_rna,
        img_pred=pred_bin,
        ct=ct,
    )
