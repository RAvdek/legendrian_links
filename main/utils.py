import json
import logging
import random
import string


with open('links.json') as f:
    LINKS = json.loads(f.read())


# Logging helpers


def get_logger(name):
    logging.basicConfig(format='%(asctime)s|%(name)s|%(levelname)s|%(message)s')
    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)
    return logger


LOG = get_logger(__name__)


def log_start_stop(func):
    """Simple function wrapper to log starting and ending times of a function

    :param func:
    :return: wrapped func
    """
    def timed(*args, **kwargs):
        LOG.info(f"Starting execution {func.__name__}")
        result = func(*args, **kwargs)
        LOG.info(f"Ending execution {func.__name__}")
        return result
    return timed


def tiny_id(k=4):
    """Generate an id string. Functional, but less safe than UUID. See https://stackoverflow.com/a/56398787.

    :param k: length of the id
    :return: str
    """
    return ''.join(random.choices(string.ascii_letters + string.digits, k=k))


# Methods for manipulating lists


def is_nonrepeating(array):
    """Are the elements in a tuple unique? Will have to do this over and over...
    Faster and less memory expensive than `len(set(tuple)) == len(tuple)` for large arrays

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
    """Get unique elements from list. Faster than list(set(array)) for large arrays

    :param array: list
    :return: list
    """
    output = list()
    for x in array:
        if x not in output:
            output.append(x)
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
