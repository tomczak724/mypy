
import numpy as np
from collections import namedtuple


def readcat(name, dtype=float, delimiter=None, comments='#'):
    '''
    Reads an ascii catalog into a named tuple. It assumes that the
    row immediately prior to the data has the names of each column.

    dtype  ------>  data type for columns
    delimiter --->  delimiter values that separates data entries
    comments  --->  character that indicates comment lines
    ''' 
    dat = open(name, 'r')

    ###  IDENTIFYING HEADER PARAMETERS
    lines = dat.readlines()
    for i, line in enumerate(lines):
        if lines[i+1][0] != comments:
            header_params = line[1:].split(delimiter)
            break
    del dat, lines

    custom_catalog = namedtuple('custom_catalog', header_params)
    catalog = custom_catalog(*np.loadtxt(name, dtype=dtype, delimiter=delimiter, unpack=1))
    return catalog



def read_sexcat(name):
    '''
    Reads an ascii catalog from sextractor into an named tuple.
    '''
    fopen = open(name, 'r')
    lines = fopen.readlines()
    fopen.close()

    i, colnames = 0, []
    while lines[i].split()[0] == '#':
        colnames.append(lines[i].split()[2].lower())
        i += 1

    sexcat = namedtuple('sexcat', colnames)
    return sexcat(*np.loadtxt(name, unpack=1))



def read_zout(name):
    '''
    #  A quick and dirty script that reads in EAZY
    #  zout files into a named tuple.
    ''' 
    dat = open(name, 'r')
    header = dat.readline()
    header_params = header[1:].split(None)
    del dat

    custom_catalog = namedtuple('custom_catalog', header_params)
    catalog = custom_catalog(*np.loadtxt(name, unpack=1))
    return catalog




