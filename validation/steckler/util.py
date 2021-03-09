import numpy as np


def print_dict(dictionary):
    for item in dictionary:
        print("key: \'{}\', value: {}".format(item, dictionary[item]))


def create_2D_array_from_dict(dict_coordinates):
    two_dimensional = []
    for v in dict_coordinates.values():
        two_dimensional.append(v)
    return two_dimensional
