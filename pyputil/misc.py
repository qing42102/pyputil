import copy
import json
import os
import time
from dataclasses import field
from types import SimpleNamespace

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


def accept_dict_args(func):
    def wrapper(args):
        if type(args) == dict:
            return func(SimpleNamespace(**args))
        else:
            return func(args)

    return wrapper
