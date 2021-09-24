import glob
import os
import random as rn
from pathlib import Path

import numpy as np
import tensorflow as tf
from PIL import Image
from keras import backend as K
from keras import regularizers
from keras.layers import Dropout, Conv2D, Input, Activation, BatchNormalization, add
from keras.models import Model
from skimage.io import imread
from tensorflow.keras.layers import Layer
from tqdm import tqdm

from genegram.shared import ROOT


# Data processing functions

# set each gray pixel of network prediction to black/white
# according to threshold coeff
def binarize_output(img, coeff=0.6):
    size = len(img)
    for i in range(size):
        for j in range(size):
            if i != j:
                if img[i][j] > 255 * coeff:
                    img[i][j] = 255
                else:
                    img[i][j] = 0
    return img


# load one image in format accepted by network (required for tests)
def load_image(f):
    img_arr = []
    img = imread(f, as_gray=True)
    n = img.shape[0]
    img_arr.append(img)
    return np.array(img_arr).reshape(1, n, n, 1)


# get model output on input file f
def predict_img(f, model):
    inp = load_image(f)
    res = model.predict(inp)
    k = len(inp[0])
    img = np.array(Image.new("L", (k, k), (255)))
    for i in range(k):
        for j in range(k):
            img[i][j] = min(res[0][i][j][0], 255)
    return img


def img2ct(img, meta, seq):
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


def predict(input_dir: Path, output_dir: Path, weights: Path):
    # setup environment
    os.environ["PYTHONHASHSEED"] = "0"
    np.random.seed(42)
    rn.seed(12345)
    session_conf = tf.compat.v1.ConfigProto(
        intra_op_parallelism_threads=1, inter_op_parallelism_threads=1
    )
    tf.random.set_seed(1234)
    sess = tf.compat.v1.Session(
        graph=tf.compat.v1.get_default_graph(), config=session_conf
    )
    K.set_session(sess)

    # load model and predict image for each sample in input_dir
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

    # create dir for `output_dir`
    data_dir = Path(output_dir).resolve()
    data_dir.mkdir(parents=True, exist_ok=True)

    files = glob.glob(str(input_dir) + "/*.png")
    for file in tqdm(files, desc="Predicting"):
        img = predict_img(file, model)
        img = binarize_output(img)
        # ...some other post-processing
        Image.fromarray(img, "L").save(file.replace(input_dir, output_dir))

    K.clear_session()
