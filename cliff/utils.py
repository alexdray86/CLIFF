import numpy as np
import math

def sigmoid(x):
    return 1 / (1 + math.exp(-x) )

def kl_divergence(p, q):
    return np.sum(np.where(p != 0, p * np.log(p / q), 0))

def rmse(x1,x2):
    MSE = np.square(np.subtract(x1,x2)).mean()
    return math.sqrt(MSE)
