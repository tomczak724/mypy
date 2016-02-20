
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



def readzout(name):
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



class res_filter:
    def __init__(self, header):
        self.header = header
        self.data = []
        self.lams = []
        self.trans = []


def read_res(res):
    '''
    #  Outputs the filter transmission curves from
    #  a *.res file commonly used with EAZY/FAST
    #
    #  res --- full path name to the *.res file
    #
    #  EXAMPLE:
    #    >>> res = '../Filters/FILTER.RES.v6.R300'
    #    >>> mypy.sps.read_filter_from_res(10, res)
    #    array([[ 2.31000e+03   0.00000e+00],
    #           [ 2.31500e+03   0.00000e+00],
    #              ...
    #           [ 9.34500e+03   0.00000e+00],
    #           [ 9.35000e+03   0.00000e+00]])
    '''

    dat = open(res, 'r')
    lines = dat.readlines()

    filters = []

    line0 = lines[0]
    n = int(line0.split()[0])
    filters.append(res_filter(line0))

    cnt = 0
    for line in lines[1:]:
        
        if cnt < n:
            cnt, lam, tran = line.split()
            cnt = int(cnt)
            lam = float(lam)
            tran = float(tran)
            filters[-1].data.append([lam, tran])
            filters[-1].lams.append(lam)
            filters[-1].trans.append(tran)

        else:
            cnt = 0
            n = int(line.split()[0])
            filters[-1].data = np.array(filters[-1].data)
            filters[-1].lams = np.array(filters[-1].lams)
            filters[-1].trans = np.array(filters[-1].trans)
            filters.append(res_filter(line))

    filters[-1].data = np.array(filters[-1].data)
    filters[-1].lams = np.array(filters[-1].lams)
    filters[-1].trans = np.array(filters[-1].trans)
    dat.close()
    
    return filters


