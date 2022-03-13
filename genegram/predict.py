from pathlib import Path

import numpy as np
import tensorflow as tf
from keras import backend as K
from keras import regularizers
from keras.layers import Dropout, Conv2D, Input, Activation, BatchNormalization, add
from keras.models import Model
from tensorflow.keras.layers import Layer

__all__ = [
    "predict",
    "setup_model",
    "clear_session",
]


def clear_session():
    K.clear_session()


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
    return model


def predict(image: np.ndarray, model: Model) -> np.ndarray:
    """Get `model` output on input `pil_image`

    Parameters
    ----------
    image: np.ndarray
        Input image
    model: Model
        Prediction model

    Returns
    -------
    array: np.ndarray
        Prediction result
    """
    n = image.shape[0]
    inputs = np.array([image]).reshape((1, n, n, 1))

    prediction = model.predict(inputs)

    result = np.full(shape=(n, n), fill_value=255, dtype=np.uint8)
    for i in range(n):
        for j in range(n):
            if prediction[0][i][j][0] < 255:
                result[i][j] = prediction[0][i][j][0]

    return result
