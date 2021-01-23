import math as m
import numpy as np

OVERFLOW_UP_LIMIT = 1e6
OVERFLOW_DOWN_LIMIT = -1e6
OVERFLOW_DECIMAL = 1e-4
NAN_OPERATE = True


# Need to implemet all numpy type and all
# overflow_up=OVERFLOW_UP_LIMIT,overflow_down=0,overflow_decim=OVERFLOW_DECIMAL):
def overflow(scalar, overflow_up=1e6, overflow_down=-1e6, overflow_decim=1e-4):
    type_of_scalar = str(type(scalar))

    if type_of_scalar == "<class 'numpy.ndarray'>" or type_of_scalar == "<class 'array'>" or type_of_scalar == "<class 'list'>":
        if len(scalar) == 2:
            scalar[0] = overflow(scalar[0], overflow_up=overflow_up,
                                 overflow_down=overflow_down, overflow_decim=overflow_decim)
            scalar[1] = overflow(scalar[1], overflow_up=overflow_up,
                                 overflow_down=overflow_down, overflow_decim=overflow_decim)
        else:
            for i in range(len(scalar)):
                scalar[i] = overflow(scalar[i], overflow_up=overflow_up,
                                     overflow_down=overflow_down, overflow_decim=overflow_decim)
        return scalar
    elif type_of_scalar == "<class 'int'>" or type_of_scalar == "<class 'float'>" or type_of_scalar == "<class 'numpy.float64'>":
        if scalar > overflow_up:
            return overflow_up
        elif scalar < overflow_down:
            return overflow_down
        elif abs(scalar) < overflow_decim:
            return overflow_decim
        elif (m.isnan(scalar) or np.isnan(scalar)) and NAN_OPERATE:
            return overflow_up
        else:
            return scalar
    else:
        return scalar
    # return scalar


def up_overflow(scalar):
    type_of_scalar = type(scalar)
    if type_of_scalar == "'numpy.ndarray'" or type_of_scalar == "array":
        for n in scalar:
            up_overflow(n)
    elif type_of_scalar == "int":
        return scalar if scalar > OVERFLOW_UP_LIMIT else scalar
    else:
        return scalar


def has_overflown(scalar):
    return scalar < OVERFLOW_UP_LIMIT


def down_overflow(scalar):
    type_of_scalar = type(scalar)
    if type_of_scalar == "'numpy.ndarray'" or type_of_scalar == "array":
        for n in scalar:
            down_overflow(n)
    elif type_of_scalar == "int":
        return scalar if scalar < OVERFLOW_DOWN_LIMIT else scalar
    else:
        return scalar


def has_down_overflown(scalar):
    return scalar < OVERFLOW_DOWN_LIMIT
