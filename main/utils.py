import json
import logging


with open('links.json') as f:
    LINKS = json.loads(f.read())


def get_logger(name):
    logging.basicConfig(format='%(asctime)s|%(name)s|%(levelname)s|%(message)s')
    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)
    return logger


def is_nonrepeating(array):
    """Are the elements in a tuple unique? Will have to do this over and over...
    Faster and less memory expensive than `len(set(tuple)) == len(tuple)`

    :param array: list
    :return: Bool
    """
    known = []
    while len(array) > 0:
        if array[0] in known:
            return False
        else:
            known.append(array[0])
            array = array[1:]
    return True


def unique_elements(array):
    """Get unique elements from list. Faster than list(set(array))

    :param array: list
    :return: list
    """
    output = set()
    for x in array:
        output.add(x)
    return list(output)


def rotate(array, n):
    """Rotates a tuple so that the nth entry is in the 0th slot
    rotate([0, 1, 2, 3, 4], 2) = [2, 3, 4, 0, 1]

    :param array: list
    :param n: how many times to rotate
    :return: list
    """
    n %= len(array)
    return array[n:] + array[:n]
