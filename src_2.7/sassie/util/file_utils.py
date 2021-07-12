#!/usr/bin/env python
# Auther: Steven C. Howell
# Purpose: Useful tools for file and directory manipulation
# Created: 10/09/2014
# $Id: file_utils.py 3036 2016-03-01 18:58:12Z schowell $

import os
import errno

def mkdir_p(path):
    '''
    make directory recursively
    adapted from http://stackoverflow.com/questions/600268/
    '''
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


def tail(f, n=10):
    '''
    return the last n lines of f
    adapted from: http://stackoverflow.com/questions/136168
    '''
    tail_str = 'tail -n %s %s' % (str(n), f)
    stdin, stdout = os.popen2(tail_str)
    stdin.close()
    lines = stdout.readlines()
    stdout.close()
    return lines[:]


class cd:
    """
    Context manager for changing the current working directory
    http://stackoverflow.com/questions/431684
    """

    def __init__(self, newPath):
        self.newPath = newPath

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)

