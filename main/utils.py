from contextlib import contextmanager
import json
import logging
import os
import time
import numpy as np
import random
import signal
import string
import threading
import sympy


with open('links.json') as f:
    LINKS = json.loads(f.read())


# Logging helpers


def get_logger(name):
    logging.basicConfig(format='%(asctime)s|%(name)s:%(lineno)d|%(levelname)s|%(message)s')
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


# Timeout context from https://www.jujens.eu/posts/en/2018/Jun/02/python-timeout-function/#:~:text=You%20can%20use%20signals%20and,alarm%20signal%20for%20the%20timeout.
# It only works in main thread as `signal.signal` can only make calls from the main thread.

def in_main_thread():
    return threading.current_thread() == threading.main_thread()

def os_is_posix():
    return os.name == 'posix'

@contextmanager
def timeout_ctx(secs):
    if not os_is_posix():
        yield
    if not in_main_thread():
        raise RuntimeError("Trying to use timeout_ctx outside of main thread.")
    # Register a function to raise a TimeoutError on the signal.
    signal.signal(signal.SIGALRM, raise_timeout_error)
    # Schedule the signal to be sent after ``time``.
    # Using signal.alarm(time) requires time to be int while signal.setitimer allows float.
    # This is essential as waiting 1 second can be extremely expensive.
    signal.setitimer(signal.ITIMER_REAL, secs)
    try:
        yield
    except TimeoutError:
        LOG.info(f"Timeout reached at {secs}s")
        pass
    finally:
        # Unregister the signal so it won't be triggered
        # if the timeout is not reached.
        signal.signal(signal.SIGALRM, signal.SIG_IGN)


def raise_timeout_error(signum, frame):
    raise TimeoutError

@contextmanager
def timeout_manager(secs):
    start_time = time.time()
    while time.time() - start_time < secs:
        yield
    raise RuntimeError


class Timeout(object):
    def __init__(self, seconds):
        self.seconds = seconds
    def __enter__(self):
        self.die_after = time.time() + self.seconds
        return self
    def __exit__(self, type, value, traceback):
        LOG.info(f"Timeout reached at {time}s")
        pass
    @property
    def timed_out(self):
        return time.time() > self.die_after


# Methods for manipulating lists

def prod(array):
    """
    Take a product of a list of algebraic objects.
    Annoyingly, math.prod missing from some versions of python...

    :param array: list
    :return: obj
    """
    if len(array) == 0:
        return 1
    output = array.pop(0)
    while len(array) > 0:
        output *= array.pop(0)
    return output


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


# methods for numpy arrays

def one_hot_array(i, shape):
    """Return an array of size shape with zeros everywhere except having a one at index i

    :param i: Index of the 1
    :param shape: length of the array
    :return: numpy array of vector shape
    """
    output = np.zeros(shape)
    output[i] = 1
    return output


# sympy has difficulty with equality
# https://docs.sympy.org/latest/explanation/gotchas.html#double-equals-signs

def poly_equal(p1, p2):
    return sympy.expand(p1 - p2) == 0


# math helper functions

def num_inverse(n, coeff_mod):
    if coeff_mod == 0:
        return sympy.Rational(n)**(-1)
    return pow(n, -1, coeff_mod)
