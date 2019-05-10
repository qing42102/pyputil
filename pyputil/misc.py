import copy
import json
import os
import time
from dataclasses import field

import numpy as np
import scipy.linalg

START_TIME = time.time()


def default_field(obj):
    return field(default_factory=lambda: copy.deepcopy(obj))


def print_json(**kwargs):
    print(json.dumps({
        "time": time.time() - START_TIME,
        "pid": os.getpid(),
        **kwargs
    }))


def rotation_matrix(axis: np.array, theta: float):
    theta %= 2 * np.pi
    axis = np.array(axis)

    return scipy.linalg.expm(
        np.cross(np.eye(3), axis / np.linalg.norm(axis)) * theta
    )
