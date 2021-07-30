#!/usr/bin/env python

import os.path
import setuptools


def expand_data(data_to_expand):
    """
    From data structure like setup.py data_files (see http://)
     [(directory/where/to/copy/the/file, [path/to/file/in/archive/file1, ...]), ...]
    but instead of the original struct this one accept to specify a directory in elements to copy.
    This function will generate one entry for all *content* of the directory and subdirectory
    recursively, to in fine copy the tree in archive in dest on the host
    the first level of directory itself is not include (which allow to rename it)
    :param data_to_expand:
    :type  data_to_expand: list of tuple
    :return: list of tuple
    """
    def remove_prefix(prefix, path):
        prefix = os.path.normpath(prefix)
        path = os.path.normpath(path)
        to_remove = len([i for i in prefix.split(os.path.sep) if i])
        truncated = [i for i in path.split(os.path.sep) if i][to_remove:]
        truncated = os.path.sep.join(truncated)
        return truncated

    data_struct = []

    for base_dest_dir, src in data_to_expand:
        base_dest_dir = os.path.normpath(base_dest_dir)
        for one_src in src:
            if os.path.isdir(one_src):
                for path, _, files in os.walk(one_src):
                    if not files:
                        continue
                    path_2_create = remove_prefix(one_src, path)
                    data_struct.append((os.path.join(base_dest_dir, path_2_create), [os.path.join(path, f) for f in files]))
            if os.path.isfile(one_src):
                data_struct.append((base_dest_dir, [one_src]))
    return data_struct


if __name__ == "__main__":
    setuptools.setup(
        package_data={
          "integron_finder": ["data/Models/*",
                              "data/Functional_annotation/*"]
        },
        data_files=expand_data([('share/integron_finder/doc/html', ['doc/build/html']),
                                ('share/integron_finder/doc/pdf', ['doc/build/latex/IntegronFinder.pdf'])
                               ]),
        test_suite='tests.run_tests.discover',
    )