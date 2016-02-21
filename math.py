
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
    


def tomczak_sampler(pdf, N=1):
    '''
    Draws random samples from the provided PDF.

    Parameters
    ----------

    pdf : array
        2d array of the PDF. The first and second columns
        must be x and P(x) respectively. P(x) is not required
        to be normalized beforehand.
    N : int
        Number of desired samplings.

    Returns
    -------
    sample : array
        Random samplings from the provided PDF.
    '''
    
    x, px = pdf[:,0], pdf[:,1]
    px = px / np.trapz(px, x)

    u = []
    for i in range(len(x)):
        xprime = x[np.arange(i)]
        pxprime = px[np.arange(i)]
        uprime = np.trapz(pxprime, xprime)
        u.append(uprime)
    u = np.array(u)

        
    sample = []
    urand = np.random.rand(N)

    for ui in urand:

        du = abs(u - ui)
        ind = np.where(du == min(du))[0][0]
        sample.append(x[ind])

    return np.array(sample)



