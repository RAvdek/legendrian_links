import json
import logging


with open('links.json') as f:
    LINKS = json.loads(f.read())


def get_logger(name):
    logging.basicConfig(format='%(asctime)s|%(name)s|%(levelname)s|%(message)s')
    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)
    return logger


def is_nonrepeating(tuple):
    """Are the elements in a tuple unique?

    :param tuple: list
    :return: Bool
    """
    return len(set(tuple)) == len(tuple)


def rotate(tuple, n):
    """Rotates a tuple n times. For example:
    rotate([1, 2, 3, 4], 1) = [4, 1, 2, 3]

    :param tuple: list
    :param n: how many times to rotate
    :return: list
    """
    return tuple[-n:] + tuple[:-n]
