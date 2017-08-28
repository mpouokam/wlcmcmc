# -*- coding: utf-8 -*-

import numpy as np

def vdot(v1, v2):
    return v1[:, 0] * v2[:, 0] + v1[:, 1] * v2[:, 1] + v1[:, 2] * v2[:, 2]


def vnorm(input):
    return np.sqrt((input ** 2).sum(1))

def edge_lengths(v):
    return vnorm(v[1:] - v[:-1])

def vrodrot(vectors, referencepoint, axis, theta):
    vectors = vectors - referencepoint
    costh = np.cos(theta)
    sinth = np.sin(theta)
    for i in range(len(vectors)):
        vectors[i] = vectors[i] * costh + np.cross(axis, vectors[i]) * sinth + axis * (np.dot(axis, vectors[i])) * (
            1 - costh)

    vectors = vectors + referencepoint
    return vectors
