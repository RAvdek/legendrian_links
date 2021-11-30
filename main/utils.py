import logging


def get_logger(name):
    logging.basicConfig(format='%(asctime)s|%(name)s|%(levelname)s|%(message)s')
    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)
    return logger
