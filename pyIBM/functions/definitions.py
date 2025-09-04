import numpy as np

def integer_part(value):
    return int(np.floor(value))


def kronecker_delta(x, y):
    return 1 if x == y else 0


def sign(value):
    if value > 0:
        return 1
    elif value < 0:
        return -1
    else:
        return 0
