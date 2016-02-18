
import numpy


def radec_sep(ra1, dec1, ra2, dec2):
    """
    This script takes in two pairs of RA & Dec and
    returns the separation between them in arcseconds
    """
    if type(ra1) == list: ra1 = np.array(ra1)
    if type(ra2) == list: ra2 = np.array(ra2)
    if type(dec1) == list: dec1 = np.array(dec1)
    if type(dec2) == list: dec2 = np.array(dec2)

    sep = numpy.sqrt((ra1 - ra2)**2 * numpy.cos((dec1 + dec2) * numpy.pi / 360)**2 + (dec1 - dec2)**2) * 3600
    return sep


def lookback_time(z, Om=0.3, Ok=0, Olam=0.7, h=0.7):
    """
    This qwik script returns the Hubble time and
    lookback time for the input cosmo. parameters.
    Units are Gyr.

    OUTPUT = [  age(z=0) ,  age(z=z) ,  lookback  ]
    """
    th = (3.09e17) / h   # HUBBLE TIME IN SECONDS
    zrange = numpy.linspace(0, z, 10000)
    Ez = numpy.sqrt( Om*(1+zrange)**3 + Ok*(1+zrange)**2 + Olam )
    tl = th * numpy.trapz( 1.0/(1+zrange)/Ez , zrange )
    return pl.array( [ round(th/365/24/60/60/10**9,3), \
                       round(tl/365/24/60/60/10**9,3), \
                       round((th-tl)/365/24/60/60/10**9,3) ] )


def vel_recess(z):
    """
    This script takes a redshift as input and
    returns a corresponding recessional velocity
    (km/s) using the relativistic doppler shift.
    """
    return ((z+1)**2 - 1)*(2.998e5) / ((z+1)**2 + 1)



def ra(hour, minute=None, second=None, delimiter=':'):
    """
    This script attempts to convert RA between
    degrees and sexadecimal.

    Examples:

    >>> ra = '05:13:54.33'
    >>> mypy.ra(ra)
    >>>   78.476375
    >>>
    >>> mypy.ra(5, 13, 54.33)
    >>>   78.476375
    >>>
    >>> mypy.ra(78.476375)
    >>>   (5, 13, 54.33)
    """

    if type(hour) == str:
    	h, m, s = hour.split(delimiter)
    	return 15 * ( float(h) + float(m)/60. + float(s)/3600.)

    if min != None:
    	return 15 * ( hour + minute/60. + second/3600. )
    else:
        h = int(hour/15. )
        m = int((hour/15. - h) * 60)
        s = ((hour/15. - h) * 60 - m) * 60
        return h,m,s



def dec(degree, arcminute=None, arcsecond=None, delimiter=':'):
	"""
	This script attempts to convert Declination
	between degrees and sexadecimal.

	Examples:

	>>> dec = '-04:45:32.12'
	>>> mypy.dec(dec)
	>>>   -4.758922
	>>>
	>>> mypy.dec(-4, 45, 32.12)
	>>>   -4.758922
	>>>
	>>> mypy.dec(-4.758922)
	>>>   (-4, 45, 32.12)
	"""

	if type(degree) == str:
		d, amin, asec = degree.split(delimiter)
		if float(d) < 0:
			return float(d) - float(amin)/60. - float(asec)/3600.
		else:
			return float(d) + float(amin)/60. + float(asec)/3600.


	if arcminute != None:
		if degree < 0:
			return degree - arcminute/60. - arcsecond/3600.
		else:
			return degree + arcminute/60. + arcsecond/3600.
	else:
		d = int(abs(degree))
		amin = int((abs(degree) - d) * 60)
		asec = ((abs(degree) - d) * 60 - amin) * 60
		if degree < 0:
			return -d, amin, asec
		else:
			return d, amin, asec


















