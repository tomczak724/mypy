import numpy
import astropy
from matplotlib import pyplot


def synphot(spectrum, bandpass):
    """
    Returns the synthetic photometry by convolving the provided spectrum with 
    the provided filter bandpass and integrating the flux.

    Parameters
    ----------
    spectrum : array_like, shape (N, 2)
        2D array representing the spectrum to be integrated. First and second
        columns must be wavelengths and fluxes respectively.

    bandpass : array_like, shape (N, 2)
        2D array representing the bandpass transmission curve. First and 
        second columns must be wavelengths and transmission respectively.

    Returns
    -------
    flux : float
        The integrated flux of ``spectrum'' through ``bandpass''

    Notes
    -----
        synphot is short for "synthetic photometry" which is another term
        for the integrated flux of a spectrum through a bandpass.
    """

    ###  normalizing bandpass        
    bandpass_normed = bandpass.copy()
    bandpass_normed[:,1] /= numpy.trapz(bandpass[:,1], bandpass[:,0])

    ###  convolving spectrum with normalized bandpass
    spectrum_interped = numpy.interp(bandpass_normed[:,0], \
                                     spectrum[:,0], \
                                     spectrum[:,1])
    fconvolved = bandpass_normed[:,1] * spectrum_interped

    flux = numpy.trapz(fconvolved, bandpass_normed[:,0])

    return flux


def uvj_select(uv, vj, z):
    """
    Classify galaxies as Star-Forming (1) or Quiescent (0) based on the 
    UVJ selection criteria of Whitaker et al. (2011), ApJ, 735, 86.

    Parameters
    ----------
    uv, vj : float_like or array_like
        Rest-frame (U-V) and (V-J) colors of galaxies. Must have same
        dimensions.

    z : float_like or array_like
        Redshift(s) of input galaxies.

    Returns
    -------
    sfqu : float_like or array_like
        Binary classification of galaxies in UVJ space:
        1 = Star-Forming
        0 = Quiescent

    Notes
    -----
        The vertical boundary at (V-J) ~ 1.5 has been removed, making
        the selection baoundary more similar to that of the related 
        NUV-r-J classification scheme.
    """

    ###  converting inputs into 1D arrays
    uv = numpy.array(uv, ndmin=1)
    vj = numpy.array(vj, ndmin=1)
    z  = numpy.array(z, ndmin=1)


    floor     = numpy.zeros(len(z))
    slope     = 0.88
    intercept = numpy.zeros(len(z))

    floor[z <= 2.] = 1.3
    floor[z >= 2.] = 1.2

    intercept[z <= 0.5] = 0.69
    intercept[z >= 0.5] = 0.59

    sfqu = numpy.ones(len(z))
    sfqu[(uv > floor) & (uv > slope*vj + intercept)] = 0

    if len(uv) == 1:
        sfqu = sfqu[0]

    return sfqu


def uvj_select_region(z, subplot=False):
    """
    Plots the UVJ selection boundary of Whitaker et al. (2011), ApJ, 735, 86
    at the input redshift, z. 

    Parameters
    ----------
    z : float_like
        Redshift to plot UVJ boundary for.

    subplot : matplotlib.axes._subplots.AxesSubplot
        Subplot in which to plot. If not supplied it will plot within the 
        most recently declared plot window.

    Returns
    -------
    uvj_boundary : matplotlib.lines.Line2D
        The matplotlib line object representing the boundary.

    Notes
    -----
        The vertical boundary at (V-J) ~ 1.5 has been removed, making
        the selection baoundary more similar to that of the related 
        NUV-r-J classification scheme.
    """

    ###  grabbing x- and y-axis bounds of plot window
    if subplot:
        axis = subplot.axis()
        xlo, xhi, ylo, yhi = axis
    else:
        axis = pyplot.axis()
        xlo, xhi, ylo, yhi = axis

    ###  parameters of UVJ boundary lines
    if z <= 0.5:
        floor = 1.3
        slope = 0.88
        intercept = 0.69
    elif z <= 2.:
        floor = 1.3
        slope = 0.88
        intercept = 0.59
    else:
        floor = 1.2
        slope = 0.88
        intercept = 0.59

    ###  inflection point of UVJ boundary lines
    x_inflect = (floor - intercept) / slope

    ###  defining relevant plot points for UVJ boundary lines
    if x_inflect <= xlo:
        xbounds = numpy.array([xlo, xhi])
        ybounds = slope * xbounds + intercept
    elif x_inflect >= xhi:
        xbounds = numpy.array([xlo, xhi])
        ybounds = numpy.array([floor, floor])
    else:
        xbounds = numpy.array([xlo, x_inflect, xhi])
        ybounds = slope * xbounds + intercept
        ybounds[0] = floor

    ###  plotting
    if subplot:
        uvj_boundary = subplot.plot(xbounds, ybounds, color='r', lw=2)[0]
        subplot.axis(axis)
    else:
        uvj_boundary = pyplot.plot(xbounds, ybounds, color='r', lw=2)[0]
        pyplot.axis(axis)

    return uvj_boundary


def nuvrj_select(nuvr, rj):
    """
    Classify galaxies as Star-Forming (1) or Quiescent (0) based on the 
    NUV-r-J selection criteria of Ilbert et al. (2010), ApJ, 709, 644.

    Parameters
    ----------
    nuvr, rj : float_like or array_like
        Rest-frame (NUV-r) and (r-J) colors of galaxies. Must have same
        dimensions.

    Returns
    -------
    sfqu : float_like or array_like
        Binary classification of galaxies in NUV-r-J space:
        1 = Star-Forming
        0 = Quiescent
    """

    ###  converting inputs into 1D arrays
    nuvr = numpy.array(nuvr, ndmin=1)
    rj = numpy.array(rj, ndmin=1)


    floor     = 3.1
    slope     = 3.
    intercept = 1.

    sfqu = numpy.ones(len(nuvr))
    sfqu[(nuvr > floor) & (nuvr > slope*rj + intercept)] = 0

    if len(nuvr) == 1:
        sfqu = sfqu[0]

    return sfqu


def nuvrj_select_region(subplot=False):
    """
    Plots the NUV-r-J selection boundary of Ilbert et al. (2011), ApJ, 709, 644

    Parameters
    ----------
    subplot : matplotlib.axes._subplots.AxesSubplot
        Subplot in which to plot. If not supplied it will plot within the 
        most recently declared plot window.

    Returns
    -------
    nuvrj_boundary : matplotlib.lines.Line2D
        The matplotlib line object representing the boundary.
    """

    ###  grabbing x- and y-axis bounds of plot window
    if subplot:
        axis = subplot.axis()
        xlo, xhi, ylo, yhi = axis
    else:
        axis = pyplot.axis()
        xlo, xhi, ylo, yhi = axis

    ###  parameters of NUV-r-J boundary lines
    floor     = 3.1
    slope     = 3.
    intercept = 1.

    ###  inflection point of UVJ boundary lines
    x_inflect = (floor - intercept) / slope

    ###  defining relevant plot points for UVJ boundary lines
    if x_inflect <= xlo:
        xbounds = numpy.array([xlo, xhi])
        ybounds = slope * xbounds + intercept
    elif x_inflect >= xhi:
        xbounds = numpy.array([xlo, xhi])
        ybounds = numpy.array([floor, floor])
    else:
        xbounds = numpy.array([xlo, x_inflect, xhi])
        ybounds = slope * xbounds + intercept
        ybounds[0] = floor

    ###  plotting
    if subplot:
        nuvrj_boundary = subplot.plot(xbounds, ybounds, color='r', lw=2)[0]
        subplot.axis(axis)
    else:
        nuvrj_boundary = pyplot.plot(xbounds, ybounds, color='r', lw=2)[0]
        pyplot.axis(axis)

    return nuvrj_boundary



