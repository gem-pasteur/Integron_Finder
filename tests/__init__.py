import os.path
import unittest
import platform


class IntegronTest(unittest.TestCase):

    _data_dir = os.path.normpath(os.path.join(os.path.dirname(__file__), "data"))

    def find_data(self, name):
        data_path = os.path.join(self._data_dir, name)
        if os.path.exists(data_path):
            return data_path
        else:
            raise RuntimeError("data '{}' does not exists".format(name))


def which(name, flags=os.X_OK):
    """
    Search PATH for executable files with the given name.

    :param name: the name of the executable to search
    :type name: str
    :param flags: os mod the name must have, default is executable (os.X_OK).
    :type flags: os file mode R_OK|R_OK|W_OK|X_OK
    :return: the path of the executable
    :rtype: string or None
    """
    result = None
    path = os.environ.get('PATH', None)
    if path is None:
        return result
    for p in os.environ.get('PATH', '').split(os.pathsep):
        p = os.path.join(p, name)
        if platform.system() == 'Windows':
            p += '.exe'
        if os.access(p, flags):
            result = p
            break
    return result
