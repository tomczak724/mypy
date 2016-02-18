
import numpy


def nmad(array):
    """
    Returns the Normalized Median Absolute Deviation.
    """
    return 1.483 * numpy.median(abs(numpy.array(array) - numpy.median(array)))


def gauss1d(x, A, mu, sig, offset=0):
    """
    A simple Gaussian function.

    y(x) = A * exp( -(x - mu)**2 / (2 * sig**2) ) + offset
    """
    return A * numpy.exp( -(x - mu)**2 / (2. * sig**2) ) + offset


def gauss2d(lenx, leny, A, mux, muy, sigx, sigy, offset=0):
    """
    A simple 2D Gaussian function.
    (lenx, leny) are the desired integer side-lengths.

    y(x) = A * exp( -(x-mux)**2/(2*sigx**2) - (y-muy)**2/(2*sigy**2) ) + offset
    """
    xgrid, ygrid = numpy.meshgrid(range(lenx), range(leny))
    return A * numpy.exp( -(xgrid - mux)**2 / (2. * sigx**2) - (ygrid - muy)**2 / (2. * sigy**2))
    
