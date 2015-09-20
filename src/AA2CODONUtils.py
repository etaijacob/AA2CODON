import os


__author__ = 'etai'


def tail(f, n=1):
    if os.path.isfile(f):
        lines = file(f, 'r').readlines()
        return lines[-n:]
    return None

