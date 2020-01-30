"""
Convolutional neural network architecture.
Author:
Date:
"""

import numpy as np
import tensorflow as tf

from tensorflow.keras.layers import Dense, Flatten, Conv2D
from tensorflow.keras import Model

##################

class CNNmodel(Model):
    """
    A convolutional neural network; the architecture is:
    Conv -> ReLU -> Conv -> ReLU -> Dense
    Note that we only need to define the forward pass here; TensorFlow will take
    care of computing the gradients for us.
    """
    def __init__(self):
        super(CNNmodel, self).__init__()
        # TODO complete constructor
        # First conv layer: 32 filters, each 5x5
        # Second conv layer: 16 filters, each 3x3

    def call(self, x):
        # TODO complete call method
        pass

def three_layer_convnet_test():
    """Test function to make sure the dimensions are working"""

    # Create an instance of the model
    cnn_model = CNNmodel()

    # TODO try out both the options below (all zeros and random)
    # shape is: number of examples (mini-batch size), width, height, depth
    #x_np = np.zeros((64, 32, 32, 3))
    #x_np = np.random.rand(64, 32, 32, 3)

    # call the model on this input and print the result
    output = cnn_model.call(x_np)
    print(output) # TODO what shape is this? does it make sense?

    # TODO look at the model parameter shapes, do they make sense?
    for v in fc_model.trainable_variables:
        print("Variable:", v.name)
        print("Shape:", v.shape)

def main():
    # test three layer function
    three_layer_convnet_test()

if __name__ == "__main__":
    main()
