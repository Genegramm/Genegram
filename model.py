import glob
import os
import random as rn
from pathlib import Path

import numpy as np
import tensorflow as tf
from PIL import Image
from keras import backend as K
from keras import regularizers
from keras.layers import Dropout, Conv2D, Input, Activation, BatchNormalization, \
    add
from keras.models import Model
from skimage.io import imread
from tensorflow.keras.layers import Layer

os.environ["PYTHONHASHSEED"] = "0"
np.random.seed(42)
rn.seed(12345)
session_conf = tf.compat.v1.ConfigProto(
    intra_op_parallelism_threads=1, inter_op_parallelism_threads=1
)
tf.random.set_seed(1234)
sess = tf.compat.v1.Session(graph=tf.compat.v1.get_default_graph(),
                            config=session_conf)
K.set_session(sess)

src_dir = str(Path(__file__).resolve().parent) + "/"
model_dir = src_dir
input_dir = src_dir + "in/"
output_dir = src_dir + "out/"

input_shape = (None, None, 1)


# Data processing functions


def binarize_output(
        img, coeff=0.6
):  # set each gray pixel of network prediction to black/white according to threshold coeff
    size = len(img)
    for i in range(size):
        for j in range(size):
            if i != j:
                if img[i][j] > 255 * coeff:
                    img[i][j] = 255
                else:
                    img[i][j] = 0
    return img


def load_image(
        f):  # load one image in format accepted by network (required for tests)
    img_arr = []
    img = imread(f, as_gray=True)
    n = img.shape[0]
    img_arr.append(img)
    return np.array(img_arr).reshape(1, n, n, 1)


def predict_img(f, model):  # get model output on input file f
    inp = load_image(f)
    res = model.predict(inp)
    k = len(inp[0])
    img = np.array(Image.new("L", (k, k), (255)))
    for i in range(k):
        for j in range(k):
            img[i][j] = min(res[0][i][j][0], 255)
    return img


def img2ct(img, meta, seq):  # transform image to .ct format
    return


def remove_multiplets(img):  # get image without multiplets
    return


# Model definition functions


class WeightedSum(
    Layer
):  # layer that for inputs i1, i2, i3, i4 returns (w1*i1 + w2*i2 + w3*i3 + w4*i4) / 4, where wi are trainable coeffs
    def __init__(self):
        super(WeightedSum, self).__init__()
        w_init = tf.keras.initializers.Ones()
        self.w1 = tf.Variable(
            initial_value=w_init(shape=(), dtype="float32"), trainable=True
        )
        self.w2 = tf.Variable(
            initial_value=w_init(shape=(), dtype="float32"), trainable=True
        )
        self.w3 = tf.Variable(
            initial_value=w_init(shape=(), dtype="float32"), trainable=True
        )
        self.w4 = tf.Variable(
            initial_value=w_init(shape=(), dtype="float32"), trainable=True
        )

    def call(self, input):
        input1, input2, input3, input4 = input
        return (
                       tf.multiply(input1, self.w1)
                       + tf.multiply(input2, self.w2)
                       + tf.multiply(input3, self.w3)
                       + tf.multiply(input4, self.w4)
               ) / 4


def res_unit(
        inputs, filters, kernels, activ, bn=False, dr=0
):  # residual unit definition, classical structure in ML
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


def res_network(
        inputs, units_num, filters, kernels, activ="relu", bn=False, dr=0
):  # residual network definition: repeating res units + skip sonnections(add layers)
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


def parallel_res_network(
        blocks_num, units_num, filters, kernels, activ="relu", bn=False, dr=0
):  # model that combines several residual networks
    inputs = Input(shape=input_shape)
    all_outputs = []
    for i in range(blocks_num):  # construct several resnets of same shape
        model = res_network(inputs, units_num, filters, kernels, activ, bn, dr)
        all_outputs.append(model.output)

    x = WeightedSum()(all_outputs)  # combine their outputs by weighted sum

    y = res_unit(
        x, filters, kernels, activ, bn, dr
    )  # add final res unit to the end of network
    x = add([x, y])
    x = Activation(activ)(x)
    outputs = x
    model = Model(inputs=inputs, outputs=outputs)
    return model


def main():
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
    model.load_weights(model_dir + "weights_7.h5")
    pred_dir = model_dir + "predicted/"
    if not os.path.exists(pred_dir):
        os.mkdir(pred_dir)
    files = glob.glob(input_dir + "*.png")
    for f in files:
        img = predict_img(f, model)
        img = binarize_output(img)
        # ...some other post-processing
        Image.fromarray(img, "L").save(f.replace(input_dir, pred_dir))

    K.clear_session()


if __name__ == "__main__":
    main()
