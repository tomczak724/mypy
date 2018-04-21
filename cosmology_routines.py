import numpy
import astropy


def radec_sep(ra1 ,dec1 ,ra2 ,dec2):
    """
    Returns the angular separation(s) in arcseconds between the input 
    sets of Right Ascensions and Declinations.

    Parameters
    ----------
    ra1, dec1 : float_like or array_like
        First set of input coordinates. Must be in decimal form and 
        both be either float_like or array_like with same length.

    ra2, dec2 : float_like or array_like
        Second set of input coordinates. If array_like and (ra1, dec1)
        are also array_like, must have same length as (ra1, dec1).
        Must also be in decimal form.

    Returns
    -------
    separations : float_like or ndarray
        Separations in arcseconds between the input sets of coordinates. 
        If both (ra1, dec1) and (ra2, dec2) are float_like then 
        ``separations'' is also float_like. If either set of input 
        coordinates is array_like then ``separations'' is ndarray with 
        length equal to that of the input array_like set(s).

    Notes
    -----
    This function is designed specifically for the Equatorial
    coordinate system. Do not use for, e.g., Galactic coordinates.
    """
    if type(ra1) == list: ra1 = numpy.array(ra1)
    if type(ra2) == list: ra2 = numpy.array(ra2)
    if type(dec1) == list: dec1 = numpy.array(dec1)
    if type(dec2) == list: dec2 = numpy.array(dec2)

    term1 = (ra1-ra2) * numpy.cos((dec1+dec2)*numpy.pi/360)
    term2 = (dec1-dec2)
    separations = (term1**2 + term2**2)**0.5
    return separations * 3600


def ra(h, m=None, s=None):
    """
    Converts the input Right Ascension from
    degrees to HMS format, or vice versa.

    Parameters
    ----------
    h : float_like
        Either input RA in degrees or hours portion if in HMS format

    m : float_like
        Minutes portion for HMS format

    s : float_like
        Seconds portion for HMS format


    Returns
    -------
    RA : float_like or ndarray
        Converted Right Ascension quantity
    """
    if m == s == None:
        RA = 15 * (h + m/60. + s/3600.)
    else:
        h = int(h / 15.)
        m = int((h / 15.-h) * 60)
        s = ((h/15. - h) * 60 - m) * 60
        RA = (h, m, s)

    return RA


def dec(d, amin=None, asec=None):
    """
    Converts the input Declination from
    degrees to HMS format, or vice versa.

    Parameters
    ----------
    d : float_like
        Either input Dec in degrees or d portion if in HMS format

    amin : float_like
        Arcminutes portion for HMS format

    asec : float_like
        Arcseconds portion for HMS format


    Returns
    -------
    Dec : float_like or ndarray
        Converted Declination quantity
    """
    if amin == asec == None:
        if d < 0: Dec = (d - amin/60. - asec/3600.)
        if d >= 0: Dec = (d + amin/60. + asec/3600.)
    else:
        d = int(abs(d))
        amin = int((abs(d) - d) * 60)
        asec = ((abs(d) - d) * 60 - amin) * 60
        if deg < 0: Dec = (-d, amin, asec
        if deg >= 0: Dec = (d, amin, asec

    return Dec


def recessional_velocity(z):
    """
    Returns the cosmological recessional velociy at the given redshift

    Parameters
    ----------
    z : float_like or array_like
        Input redshift(s)

    Returns
    -------
    velociy : float_like or ndarray
        Recessional velociy in km/s
    """
    numerator   = ((z+1)**2 - 1) * (2.998e5) 
    denominator = ((z+1)**2 + 1)
    velociy = numerator / denominator
    return velociy











