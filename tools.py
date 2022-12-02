import os
import numpy as np

def join_path(pathelts, endWithSlash=False):
    """Join strings into a path understandable by the system

    Args:
        pathelts (list): an array with elements which constitute path
        endWithSlash (bool, optional): Does path contain a slash at the end?. Defaults to False.

    Returns:
        str: path as a single string
    """
    o = os.sep.join(pathelts)
    if endWithSlash: o += os.sep
    return o

def append_path(currentPath, newElts, endWithSlash=False):
    """Appends a list of strings to a currently valid path

    Args:
        currentPath (str): a valid (system-dependent) path. Presumably output from join_path
        newElts (list): new elements to be added
        endWithSlash (bool, optional): Does path contain a slash at the end?. Defaults to False.

    Returns:
        str: path as a single string
    """
    if not endWithSlash: return currentPath + join_path(newElts)
    else: return currentPath + join_path(newElts) + os.sep

# BASE = ["", "Users", "morgunov"]
# newElts = ["vv", "core_excitations"]
# path = join_path(BASE, endWithSlash=True)
# print(path)
# new = append_path(path, newElts, endWithSlash=True)
# print(new)

def make_combinations(data, perm = None):
    """A generator designed to create permutations of molecules. Currently unused

    Args:
        data (list): a list of molecules which we use or not
        perm (list, optional): a list with permutations. Defaults to None.

    Yields:
        list: a subset of initial "data" list
    """
    if not perm:
        perm = []
    if not data:
        if perm:
            yield perm
    else:
        yield from make_combinations(data[1:], perm +[data[0],] )
        yield from make_combinations(data[1:], perm )

def get_next_col(col, prefix = None):
    if prefix is None:
        prefix = ''
    strings = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    if not col:
        return 'A'
    if len(col) == 1:
        if col == 'Z':
            return get_next_col(prefix, prefix=None) + get_next_col('', prefix=None)
        else:
            return prefix + strings[strings.index(col)+1]
    else:
        return get_next_col(col[1:], prefix = prefix + col[0])

def rounder(dig):
    def rounder_to_dig(float):
        return format(np.around(float, dig), f'.{dig}f')
    return rounder_to_dig