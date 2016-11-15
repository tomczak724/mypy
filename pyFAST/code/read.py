
import numpy
import subprocess
from collections import namedtuple


def readparams(name, comments='#'):

    ###  gunzip if necessary
    gzipped = 0
    if name[-3:] == '.gz':
        gzipped = 1
        subprocess.call('gunzip %s' % name, shell=1)
        name = name[:-3]

    initial_read = numpy.loadtxt(name, dtype=str, comments='#')
    for i in range(len(initial_read)):
        initial_read[i][2] = initial_read[i][2].strip()
        initial_read[i][2] = initial_read[i][2].strip("'")

    custom_catalog = namedtuple('custom_catalog', initial_read[:,0])
    catalog = custom_catalog(*initial_read[:,2])

    ###  re-gzip if necessary
    if gzipped:
        subprocess.call('gzip %s' % name, shell=1)

    return catalog

class res_filter:
    def __init__(self, header):
        self.header = header
        self.data = []
        self.lams = []
        self.trans = []
        self.l0 = 0.
def readres(res):
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
            filters[-1].data = numpy.array(filters[-1].data)
            filters[-1].lams = numpy.array(filters[-1].lams)
            filters[-1].trans = numpy.array(filters[-1].trans)
            filters.append(res_filter(line))

    filters[-1].data = numpy.array(filters[-1].data)
    filters[-1].lams = numpy.array(filters[-1].lams)
    filters[-1].trans = numpy.array(filters[-1].trans)
    dat.close()

    ###  estimating central wavelength
    for f in filters:
        areas = []
        for i in range(len(f.lams)):
            areas.append(numpy.trapz(f.trans[:i], f.lams[:i]))
        f.l0 = numpy.interp(0.5, numpy.array(areas)/max(areas), f.lams)

    
    return filters


