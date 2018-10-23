from collections import namedtuple
from pathlib import Path
import pkgutil
from importlib import import_module

"""Sequence description with fields: id strand start stop"""
SeqDesc = namedtuple('SeqDesc', ('id', 'strand', 'start', 'stop'))


def make_register():
    """

    :return:
    """
    parsers = {}

    def register(func):
        """

        :param func:
        :return:
        """
        parsers[func.__name__] = func
        return func

    def load_parsers():
        """

        :return:
        """
        for (_, module_name, _) in pkgutil.iter_modules([Path(__file__).parent]):
            import_module('.' + module_name, package=__name__)
        return parsers.copy()

    return register, load_parsers


register, load_parsers = make_register()